// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import Paper from '@material-ui/core/Paper';
import Box from '@material-ui/core/Box';
import FormGroup from '@material-ui/core/FormGroup';
import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import CircularProgress from '@material-ui/core/CircularProgress';
import { green } from '@material-ui/core/colors';
import { Formik, Form } from 'formik';
import * as Yup from 'yup';
import Backdrop from '@material-ui/core/Backdrop';
import { InputLabel } from '@material-ui/core';
import Common from './Common';
import * as Api from '../../api';
import { JOBS } from '../../constants/routes.json';
import SelectField from '../Form/SelectField';
import TextField from '../Form/TextField';
import Wizard from '../UI/Wizard';
import FileSelector from '../UI/FileSelector2';
import type { File } from '../UI/FileSelector2';
import UploadProgress from '../UI/UploadProgress';
import SwitchField from '../Form/SwitchField';
import type { CircRNAAnalysisConfig } from '../../types/analysis';
import type {
  Capabilities,
  ResponseType,
  SimpleMapType
} from '../../types/common';
import type { Job } from '../../types/jobs';
import ValidationError from '../../errors/ValidationError';
import SamplesField, {
  prepareFileArray,
  prepareSamplesArray
} from '../UI/SamplesField';
import CustomArgumentsField, {
  ProcessCustomArguments
} from '../Form/CustomArgumentsField';

type Props = {
  refreshJobs: () => void,
  redirect: mixed => void,
  pushNotification: (
    string,
    ?('success' | 'warning' | 'error' | 'info')
  ) => void,
  classes: {
    root: *,
    formControl: *,
    buttonWrapper: *,
    buttonProgress: *,
    backButton: *,
    instructions: *,
    backdrop: *
  }
};

const style = theme => ({
  root: {
    padding: theme.spacing(3, 2)
  },
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120
  },
  buttonWrapper: {
    margin: theme.spacing(1),
    position: 'relative'
  },
  buttonProgress: {
    color: green[500],
    position: 'absolute',
    top: '50%',
    left: '50%',
    marginTop: -12,
    marginLeft: -12
  },
  backButton: {
    marginRight: theme.spacing(1)
  },
  instructions: {
    marginTop: theme.spacing(1),
    marginBottom: theme.spacing(1)
  },
  backdrop: {
    zIndex: theme.zIndex.drawer + 1,
    color: '#fff'
  }
});

type State = {
  isLoading: boolean,
  fetched: boolean,
  isSaving: boolean,
  capabilities: ?Capabilities,
  genomes: SimpleMapType<string>,
  annotations: SimpleMapType<string>,
  bedAnnotations: SimpleMapType<string>,
  hasValidationErrors: boolean,
  validationErrors: *,
  isUploading: boolean,
  uploadFile: string,
  uploadedBytes: number,
  uploadedPercent: number,
  uploadTotal: number
};

class CircRNA extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.state = {
      isLoading: false,
      fetched: false,
      isSaving: false,
      capabilities: undefined,
      ...Api.Upload.ui.initUploadState(),
      hasValidationErrors: false,
      validationErrors: {},
      genomes: {},
      annotations: {},
      bedAnnotations: {}
    };
  }

  componentDidMount(): void {
    const { pushNotification } = this.props;
    const { fetched } = this.state;
    if (!fetched) {
      this.fetchData().catch(e => {
        pushNotification(`An error occurred: ${e.message}!`, 'error');
      });
    }
  }

  async fetchData(): Promise<void> {
    const { pushNotification } = this.props;
    try {
      this.setState({ isLoading: true });
      const genomes = await Api.References.fetchAllByAlgorithm('bwa');
      const annotations = await Api.Annotations.fetchAllByType('gtf');
      const bedAnnotations = await Api.Annotations.fetchAllByType('bed');
      const capabilities = await Api.Utils.refreshCapabilities();
      this.setState({
        isLoading: false,
        fetched: true,
        capabilities,
        genomes,
        annotations,
        bedAnnotations
      });
    } catch (e) {
      pushNotification(`An error occurred: ${e.message}!`, 'error');
      this.setState({ isLoading: false });
    }
  }

  getValidationSchema = () => {
    const { capabilities } = this.state;
    const cores = capabilities ? capabilities.availableCores : 1;
    return Yup.object().shape({
      code: Yup.string()
        .required()
        .matches(/^[A-Za-z0-9\-_]+$/, {
          message: 'The field must contain only letters, numbers, and dashes.'
        }),
      name: Yup.string().required(),
      paired: Yup.boolean().required(),
      inputType: Yup.string().oneOf(
        Object.keys(Api.Utils.supportedAnalysisFileTypes())
      ),
      convertBam: Yup.boolean().notRequired(),
      trimGalore: Yup.object()
        .notRequired()
        .shape({
          enable: Yup.boolean().notRequired(),
          quality: Yup.number()
            .notRequired()
            .min(1),
          length: Yup.number()
            .notRequired()
            .min(40),
          hardTrim: Yup.boolean().notRequired()
        }),
      genome: Yup.string().notRequired(),
      annotation: Yup.string().notRequired(),
      bedAnnotation: Yup.string().notRequired(),
      useFastqPair: Yup.boolean().notRequired(),
      ciriSpanningDistance: Yup.number().notRequired(),
      ciriUseVersion1: Yup.boolean().notRequired(),
      threads: Yup.number()
        .required()
        .min(1)
        .max(cores)
    });
  };

  getSteps = () => [
    'Choose type',
    'Set pipeline preferences',
    'Select references',
    'Upload files'
  ];

  getStep0 = () => {
    const { classes } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can choose sample code, analysis name, input file type, and
          sequencing strategy (single or paired-end). Sample code is used for
          further analysis, such as Differential Expression Analysis. It
          identifies the sample during the analysis therefore it should be a
          string without any spaces (only letters, numbers, and dashes).
        </Typography>
        <TextField label="Sample Code" name="code" required />
        <TextField label="Analysis Name" name="name" required />
        <SelectField
          label="Input Type"
          name="inputType"
          options={Api.Utils.supportedAnalysisFileTypes()}
          required
        />
        <SwitchField label="Are reads paired-end?" name="paired" />
      </>
    );
  };

  threadsText = values => {
    const { capabilities } = this.state;
    const allCores = capabilities ? capabilities.numCores : 1;
    const { threads } = values;
    const maxMultiple = Math.floor(allCores / 3);
    const standardMessage = `Do not select more than ${maxMultiple} cores to allow for multiple concurrent analysis.`;
    if (threads <= maxMultiple) {
      return standardMessage;
    }
    return (
      <Typography color="error" component="span">
        {standardMessage}
      </Typography>
    );
  };

  getStep1 = values => {
    const { classes } = this.props;
    const {
      inputType,
      convertBam,
      paired,
      ciriQuant,
      trimGalore: { enable },
      ciriUseVersion1
    } = values;
    const samOrBam = inputType === 'bam' || inputType === 'sam';
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can choose which steps will be included in the analysis:
          CIRIquant analysis (for paired-end sequencing), trimming, BAM to FASTQ
          conversion, and CIRI version.
        </Typography>
        {paired && <SwitchField label="Use CIRIquant?" name="ciriQuant" />}
        {!ciriQuant && samOrBam && (
          <SwitchField label="Convert BAM/SAM to FASTQ?" name="convertBam" />
        )}
        {(inputType === 'fastq' || convertBam || (ciriQuant && samOrBam)) && (
          <>
            <SwitchField label="Enable Trimming?" name="trimGalore.enable" />
            {enable && (
              <>
                <FormGroup row className={classes.formControl}>
                  <Grid
                    container
                    justify="center"
                    alignItems="flex-start"
                    spacing={1}
                  >
                    <Grid item xs={6}>
                      <TextField
                        label="Min PHREAD quality"
                        name="trimGalore.quality"
                        type="number"
                      />
                    </Grid>
                    <Grid item xs={6}>
                      <TextField
                        label="Min reads length"
                        name="trimGalore.length"
                        type="number"
                      />
                    </Grid>
                  </Grid>
                </FormGroup>
                <CustomArgumentsField
                  name="trimGalore.custom_arguments"
                  labelEnable="Enable Custom Arguments for TrimGalore?"
                />
              </>
            )}
          </>
        )}
        {!ciriQuant && (
          <CustomArgumentsField
            name="alignment_custom_arguments"
            labelEnable="Enable Custom Arguments for the BWA algorithm?"
          />
        )}
        {!ciriQuant && (
          <SwitchField label="Use CIRI v. 1?" name="ciriUseVersion1" />
        )}
        {!ciriQuant && ciriUseVersion1 && (
          <SwitchField label="Hard Trim reads?" name="trimGalore.hardTrim" />
        )}
        {!ciriQuant && (
          <SwitchField
            label="Use fastq-pair to reduce memory usage?"
            name="useFastqPair"
          />
        )}
        <TextField
          label="Number of threads"
          name="threads"
          type="number"
          helperText={this.threadsText(values)}
          required
        />
      </>
    );
  };

  getStep2 = () => {
    const { classes } = this.props;
    const { genomes, annotations, bedAnnotations } = this.state;
    return (
      <>
        <Typography className={classes.instructions}>
          Choose reference genome, and genome annotations if required.
        </Typography>
        <SelectField
          label="Reference Genome"
          name="genome"
          options={genomes}
          required
        />
        <SelectField
          label="Genome Annotation"
          name="annotation"
          options={annotations}
          required
        />
        <SelectField
          label="BED Back-Spliced Junction Annotation"
          name="bedAnnotation"
          options={bedAnnotations}
          required
        />
      </>
    );
  };

  getStep3 = values => {
    const { classes } = this.props;
    const { samples, code, inputType, paired } = values;
    const {
      hasValidationErrors,
      validationErrors,
      isUploading,
      uploadFile,
      uploadedBytes,
      uploadedPercent,
      uploadTotal
    } = this.state;
    const single = samples.length <= 1;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can add samples to the analysis and upload their files. For
          each sample, you will be also able to input a custom Sample Code. If
          you are uploading multiple samples for a batch analysis, you can also
          select a sample description file (in TSV format) that can be used for
          differential expression analysis. To start the upload process and the
          analysis, click on the &quot;Start Analysis&quot; button.
        </Typography>
        <SamplesField
          name="samples"
          code={code}
          inputType={inputType}
          paired={paired && inputType === 'fastq'}
          disabled={isUploading}
        />
        {!single && (
          <FormGroup row className={classes.formControl}>
            <Grid container alignItems="center" spacing={3}>
              <Grid item xs="auto">
                <InputLabel>Select an optional description file</InputLabel>
              </Grid>
              <Grid item xs>
                <FileSelector
                  name="descriptionFile"
                  filters={[
                    {
                      name: 'TSV files',
                      extensions: ['tab', 'tsv', 'txt']
                    }
                  ]}
                  disabled={isUploading}
                />
              </Grid>
            </Grid>
          </FormGroup>
        )}
        {hasValidationErrors && (
          <FormGroup row className={classes.formControl}>
            <Grid container alignItems="center" spacing={3}>
              <Grid item xs="auto">
                <Typography color="error" paragraph variant="caption">
                  Error log:
                </Typography>
                <pre>{JSON.stringify(validationErrors)}</pre>
              </Grid>
            </Grid>
          </FormGroup>
        )}
        <UploadProgress
          isUploading={isUploading}
          uploadFile={uploadFile}
          uploadedBytes={uploadedBytes}
          uploadedPercent={uploadedPercent}
          uploadTotal={uploadTotal}
        />
      </>
    );
  };

  getSubmitButton = () => {
    const { classes } = this.props;
    const { isSaving } = this.state;
    return (
      <>
        <Button
          type="submit"
          variant="contained"
          color="primary"
          disabled={isSaving}
        >
          Start Analysis
        </Button>
        {isSaving && (
          <CircularProgress size={24} className={classes.buttonProgress} />
        )}
      </>
    );
  };

  setSaving = (isSaving, validationErrors = {}) => {
    const hasValidationErrors = !(
      Object.keys(validationErrors).length === 0 &&
      validationErrors.constructor === Object
    );
    this.setState({
      isSaving,
      hasValidationErrors,
      validationErrors
    });
  };

  checkParameters = (
    parameters: CircRNAAnalysisConfig
  ): CircRNAAnalysisConfig => {
    const result = { ...parameters };
    if (!parameters.paired) {
      result.ciriQuant = false;
      if (!parameters.ciriUseVersion1)
        result.trimGalore = {
          ...result.trimGalore,
          hardTrim: false
        };
    } else if (parameters.ciriQuant) {
      result.useFastqPair = false;
      if (parameters.inputType === 'bam' || parameters.inputType === 'sam') {
        result.convertBam = true;
      }
    }
    return result;
  };

  createAnalysis = async (
    code: string,
    name: string,
    parameters: CircRNAAnalysisConfig,
    firstFile: File,
    secondFile: ?File
  ): Promise<Job> => {
    const data: ResponseType<Job> = await Api.Analysis.createCircRNA(
      code,
      name,
      this.checkParameters(parameters)
    );
    if (data.validationErrors) {
      this.setSaving(false, data.validationErrors);
      throw new ValidationError('Validation of input parameters failed');
    }
    if (!data.data) {
      throw new Error('Unknown error!');
    }
    const { data: job } = data;
    await Common.uploadFile(job, firstFile, this.setState.bind(this));
    await Common.uploadFile(job, secondFile, this.setState.bind(this));
    return job;
  };

  submitAnalysis = async (analysisJob: Job[], groupJob: ?Job) => {
    const { pushNotification } = this.props;
    // Submit all jobs in order of creation
    for (let i = 0, l = analysisJob.length; i < l; i += 1) {
      // eslint-disable-next-line no-await-in-loop
      await Api.Jobs.submitJob(analysisJob[i].id);
    }
    if (groupJob) {
      await Api.Jobs.submitJob(groupJob.id);
    }
    pushNotification('Analysis jobs queued!');
  };

  formSubmit = async values => {
    const { paired, inputType } = values;
    const {
      code,
      name,
      samples: samplesData,
      descriptionFile: description,
      trimGalore: trimGaloreParameters,
      alignment_custom_arguments: alignmentCustomArguments,
      ...tempParams
    } = values;
    const {
      custom_arguments: trimGaloreCustomArguments,
      ...trimGalore
    } = trimGaloreParameters;
    const params = ProcessCustomArguments(
      {
        ...tempParams,
        trimGalore: ProcessCustomArguments(
          trimGalore,
          trimGaloreCustomArguments,
          'custom_arguments'
        )
      },
      alignmentCustomArguments,
      'alignment_custom_arguments'
    );
    const isPairedFiles = paired && inputType === 'fastq';
    const { pushNotification, redirect, refreshJobs } = this.props;
    const files = prepareFileArray(samplesData);
    const validLength = Common.checkLength(
      files,
      isPairedFiles,
      pushNotification
    );
    if (validLength < 0) return;
    const single = validLength === 1;

    const samples = prepareSamplesArray(samplesData);

    const descriptionFile =
      description && description[0] ? description[0] : null;

    this.setSaving(true);
    try {
      const jobs = [];
      for (let i = 0; i < samples.length; i += 1) {
        const sample = samples[i];
        const [firstFile, secondFile] = files[i];
        if (firstFile && (!isPairedFiles || (isPairedFiles && secondFile))) {
          const idx = i + 1;
          const sampleCode = sample.code
            ? sample.code
            : `${code}${single ? '' : `_${idx}`}`;
          const sampleName = `${name}${single ? '' : ` - Sample ${idx}`}`;
          pushNotification(`Creating job ${sampleName}!`, 'info');
          jobs.push(
            // eslint-disable-next-line no-await-in-loop
            await this.createAnalysis(
              sampleCode,
              sampleName,
              {
                ...params,
                firstInputFile: firstFile.name,
                secondInputFile:
                  isPairedFiles && secondFile ? secondFile.name : null
              },
              firstFile,
              isPairedFiles ? secondFile : null
            )
          );
          pushNotification(`Job ${sampleName} created!`, 'success');
        }
      }
      if (!single) pushNotification('Analysis jobs created!');
      const groupJob = await Common.createGroup(
        single,
        code,
        name,
        jobs,
        descriptionFile,
        pushNotification,
        this.setSaving,
        this.setState.bind(this)
      );
      if (groupJob !== undefined) {
        await this.submitAnalysis(jobs, groupJob);
        refreshJobs();
        redirect(JOBS);
      }
    } catch (e) {
      pushNotification(`An error occurred: ${e.message}`, 'error');
      if (!(e instanceof ValidationError)) {
        this.setSaving(false);
      }
    }
  };

  render() {
    const { classes } = this.props;
    const { validationErrors, isLoading } = this.state;
    const steps = this.getSteps();
    return (
      <>
        <Box>
          <Paper className={classes.root}>
            <Typography variant="h5" component="h3">
              CircRNA Analysis
            </Typography>
            <Formik
              initialValues={{
                code: '',
                name: '',
                paired: false,
                inputType: 'fastq',
                convertBam: false,
                trimGalore: {
                  enable: true,
                  quality: 20,
                  length: 40,
                  hardTrim: false,
                  custom_arguments: {
                    enable: false,
                    value: ''
                  }
                },
                alignment_custom_arguments: {
                  enable: false,
                  value: ''
                },
                genome: 'Human_hg19_genome',
                annotation: 'Human_hg19_gencode_19_gtf',
                bedAnnotation: 'Human_hg19_circRNAs_bed',
                threads: 1,
                useFastqPair: false,
                ciriSpanningDistance: 200000,
                ciriUseVersion1: false,
                samples: [
                  {
                    code: '',
                    first: undefined,
                    second: undefined
                  }
                ],
                descriptionFile: undefined
              }}
              initialErrors={validationErrors}
              validationSchema={this.getValidationSchema()}
              onSubmit={v => {
                this.formSubmit(v).catch(() => false);
              }}
            >
              {({ values }) => (
                <Form>
                  <Wizard steps={steps} submitButton={this.getSubmitButton}>
                    <div>{this.getStep0()}</div>
                    <div>{this.getStep1(values)}</div>
                    <div>{this.getStep2()}</div>
                    <div>{this.getStep3(values)}</div>
                  </Wizard>
                </Form>
              )}
            </Formik>
          </Paper>
        </Box>
        <Backdrop className={classes.backdrop} open={isLoading}>
          <CircularProgress color="inherit" />
        </Backdrop>
      </>
    );
  }
}

export default withStyles(style)(CircRNA);
