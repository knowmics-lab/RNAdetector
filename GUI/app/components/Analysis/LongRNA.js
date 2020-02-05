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
import { Formik, Form, FieldArray } from 'formik';
import * as Yup from 'yup';
import Backdrop from '@material-ui/core/Backdrop';
import { InputLabel } from '@material-ui/core';
import IconButton from '@material-ui/core/IconButton';
import Icon from '@material-ui/core/Icon';
import * as Api from '../../api';
import { JOBS } from '../../constants/routes';
import SelectField from '../Form/SelectField';
import TextField from '../Form/TextField';
import Wizard from '../UI/Wizard';
import FileSelector from '../UI/FileSelector';
import type { File } from '../UI/FileSelector';
import UploadProgress from '../UI/UploadProgress';
import SwitchField from '../Form/SwitchField';
import type { LongRNAAnalysisConfig } from '../../types/analysis';
import type { SimpleMapType } from '../../types/common';
import type { Job } from '../../types/jobs';

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

const ALGORITHMS = {
  salmon: 'Salmon',
  tophat: 'Tophat',
  hisat2: 'Hisat 2'
};

const COUNTING_ALGORITHMS = {
  htseq: 'HT-seq',
  'feature-counts': 'Feature Counts',
  salmon: 'Salmon'
};

type State = {
  isLoading: boolean,
  fetched: boolean,
  isSaving: boolean,
  files: (?File)[][],
  descriptionFile: ?File,
  genomes: SimpleMapType<SimpleMapType<string>>,
  annotations: SimpleMapType<string>,
  validationErrors: *,
  isUploading: boolean,
  uploadFile: string,
  uploadedBytes: number,
  uploadedPercent: number,
  uploadTotal: number
};

class LongRNA extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.state = {
      isLoading: false,
      fetched: false,
      isSaving: false,
      ...Api.Upload.ui.initUploadState(),
      files: [[null, null]],
      descriptionFile: null,
      validationErrors: {},
      genomes: {},
      annotations: {}
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
      const algoKeys = Object.keys(ALGORITHMS);
      const algoValues = await Promise.all(
        algoKeys.map(a =>
          Api.References.fetchAllByAlgorithm(a === 'hisat2' ? 'hisat' : a)
        )
      );
      const annotations = await Api.Annotations.fetchAllByType('gtf');
      this.setState({
        isLoading: false,
        fetched: true,
        genomes: Object.fromEntries(algoKeys.map((a, i) => [a, algoValues[i]])),
        annotations
      });
    } catch (e) {
      pushNotification(`An error occurred: ${e.message}!`, 'error');
      this.setState({ isLoading: false });
    }
  }

  getValidationSchema = () =>
    Yup.object().shape({
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
            .min(1)
        }),
      algorithm: Yup.string()
        .notRequired()
        .oneOf(Object.keys(ALGORITHMS)),
      countingAlgorithm: Yup.string()
        .notRequired()
        .oneOf(Object.keys(COUNTING_ALGORITHMS)),
      genome: Yup.string().notRequired(),
      transcriptome: Yup.string().notRequired(),
      annotation: Yup.string().notRequired(),
      threads: Yup.number()
        .required()
        .min(1)
        .max(Api.Utils.cpuCount())
    });

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
        <TextField
          label="Number of threads"
          name="threads"
          type="number"
          helperText={`Do not select more than ${Math.floor(
            Api.Utils.cpuCount() / 3
          )} cores if you wish to submit multiple analysis.`}
          required
        />
      </>
    );
  };

  getStep1 = values => {
    const { classes } = this.props;
    const {
      inputType,
      algorithm,
      convertBam,
      trimGalore: { enable }
    } = values;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can choose which steps will be included in the analysis:
          trimming, BAM to FASTQ conversion, alignment and counting, or
          quantification.
        </Typography>
        {(inputType === 'bam' || inputType === 'sam') && (
          <SwitchField label="Convert BAM/SAM to FASTQ?" name="convertBam" />
        )}
        {(inputType === 'fastq' || convertBam) && (
          <>
            <SwitchField label="Enable Trimming?" name="trimGalore.enable" />
            {enable && (
              <FormGroup row className={classes.formControl}>
                <Grid
                  container
                  justify="center"
                  alignItems="center"
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
            )}
          </>
        )}
        <SelectField
          label="Aligment/Quantification Algorithm"
          name="algorithm"
          options={ALGORITHMS}
        />
        {(algorithm === 'tophat' || algorithm === 'hisat2') && (
          <SelectField
            label="Counting Algorithm"
            name="countingAlgorithm"
            options={COUNTING_ALGORITHMS}
          />
        )}
      </>
    );
  };

  getStep2 = values => {
    const { classes } = this.props;
    const { genomes, annotations } = this.state;
    const { algorithm, countingAlgorithm } = values;
    return (
      <>
        <Typography className={classes.instructions}>
          Choose reference genome/transcriptome, and genome annotations if
          required.
        </Typography>
        {(algorithm === 'tophat' || algorithm === 'hisat2') && (
          <SelectField
            label="Reference Genome"
            name="genome"
            options={genomes[algorithm]}
            required
          />
        )}
        {(algorithm === 'salmon' || countingAlgorithm === 'salmon') && (
          <SelectField
            label="Reference Transcriptome"
            name="transcriptome"
            options={genomes.salmon}
            required
          />
        )}
        {((algorithm !== 'salmon' && countingAlgorithm !== 'salmon') ||
          algorithm === 'tophat') && (
          <SelectField
            label="Genome Annotation"
            name="annotation"
            options={annotations}
            required
          />
        )}
      </>
    );
  };

  updateFile = (state: State, i: number, j: number, f: ?File) => {
    const files = [...state.files];
    const fileI = [...files[i]];
    fileI[j] = f;
    files[i] = fileI;
    return {
      files
    };
  };

  makeHandleFileAdd = (j: number, i: number) => f =>
    this.setState(oldState => this.updateFile(oldState, i, j, f[0]));

  makeHandleFileRemove = (j: number, i: number) => () =>
    this.setState(oldState => this.updateFile(oldState, i, j, null));

  handleDescriptionFileAdd = (descriptionFile: File[]) =>
    this.setState({
      descriptionFile: descriptionFile[0]
    });

  handleDescriptionFileRemove = () =>
    this.setState({
      descriptionFile: null
    });

  addHandle = helpers => e => {
    e.preventDefault();
    helpers.push({ code: '' });
    this.setState(oldState => ({
      files: [...oldState.files, [null, null]]
    }));
  };

  removeHandle = (helpers, i) => e => {
    e.preventDefault();
    helpers.remove(i);
    this.setState(oldState => {
      const files = [...oldState.files];
      files.splice(i, 1);
      return {
        files
      };
    });
  };

  makeSampleForm = (i, single, values, helpers) => {
    const { classes } = this.props;
    const { code, inputType, paired } = values;
    const { isUploading } = this.state;
    return (
      <FormGroup row className={classes.formControl} key={`sample-${i}`}>
        <Grid container justify="space-around" alignItems="center" spacing={3}>
          {!single && (
            <Grid item xs>
              <TextField
                label="Sample Code"
                name={`samples.${i}.code`}
                placeholder={`${code}_${i + 1}`}
              />
            </Grid>
          )}
          <Grid item xs>
            <FileSelector
              onFileRemove={this.makeHandleFileRemove(0, i)}
              onFileAdd={this.makeHandleFileAdd(0, i)}
              filters={Api.Utils.analysisFileExtensions(inputType)}
              disabled={isUploading}
            />
          </Grid>
          {paired && (
            <Grid item xs>
              <FileSelector
                onFileRemove={this.makeHandleFileRemove(1, i)}
                onFileAdd={this.makeHandleFileAdd(1, i)}
                filters={Api.Utils.analysisFileExtensions(inputType)}
                disabled={isUploading}
              />
            </Grid>
          )}
          {!single && (
            <Grid item xs={1}>
              <IconButton onClick={this.removeHandle(helpers, i)}>
                <Icon className="fas fa-trash" />
              </IconButton>
            </Grid>
          )}
        </Grid>
      </FormGroup>
    );
  };

  getStep3 = values => {
    const { classes } = this.props;
    const { samples } = values;
    const {
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
        <FieldArray
          name="samples"
          render={helpers => (
            <>
              {samples.map((s, i) =>
                this.makeSampleForm(i, single, values, helpers)
              )}
              <FormGroup row className={classes.formControl}>
                <Grid
                  container
                  direction="row-reverse"
                  alignItems="center"
                  spacing={3}
                >
                  <Grid item xs="auto">
                    <Button
                      variant="outlined"
                      onClick={this.addHandle(helpers)}
                    >
                      <Icon className="fas fa-plus" /> Add another sample
                    </Button>
                  </Grid>
                </Grid>
              </FormGroup>
            </>
          )}
        />
        {!single && (
          <FormGroup row className={classes.formControl}>
            <Grid container alignItems="center" spacing={3}>
              <Grid item xs="auto">
                <InputLabel>Select an optional description file</InputLabel>
              </Grid>
              <Grid item xs>
                <FileSelector
                  onFileRemove={this.handleDescriptionFileRemove}
                  onFileAdd={this.handleDescriptionFileAdd}
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
    this.setState({
      isSaving,
      validationErrors
    });
  };

  uploadFile = async (job: Job, file: File) => {
    Api.Upload.ui.uploadStart(this.setState.bind(this), file.name);
    await Api.Upload.upload(
      job,
      file.path,
      file.name,
      file.type,
      Api.Upload.ui.makeOnProgress(this.setState.bind(this))
    );
    Api.Upload.ui.uploadEnd(this.setState.bind(this));
  };

  createAnalysis = async (
    code: string,
    name: string,
    parameters: LongRNAAnalysisConfig,
    firstFile: File,
    secondFile: ?File
  ): Promise<Job> => {
    const data = await Api.Analysis.createLongRNA(code, name, parameters);
    if (data.validationErrors) {
      this.setSaving(false, data.validationErrors);
      throw new Error('Validation of input parameters failed');
    }
    const { data: job } = data;
    await this.uploadFile(job, firstFile);
    if (secondFile) await this.uploadFile(job, secondFile);
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

  createGroup = async (
    single: boolean,
    code: string,
    name: string,
    jobs: Job[],
    descriptionFile: ?File
  ): Promise<?Job> => {
    if (single) return null;
    const { pushNotification } = this.props;
    const data = await Api.Analysis.createSampleGroup(
      code,
      name,
      jobs,
      descriptionFile ? descriptionFile.name : undefined
    );
    if (data.validationErrors) {
      pushNotification(
        'Errors occurred during validation of input parameters. Please review the form!',
        'warning'
      );
      this.setSaving(false, data.validationErrors);
      return null;
    }
    const { data: job } = data;
    if (descriptionFile) {
      await this.uploadFile(job, descriptionFile);
    }
    pushNotification('Sample group created!');
    return job;
  };

  formSubmit = async values => {
    const { paired } = values;
    const { code, name, samples, ...params } = values;
    const { pushNotification, redirect, refreshJobs } = this.props;
    const { files, descriptionFile } = this.state;
    const validLength = files.filter(
      f => f[0] !== null && (!paired || (paired && f[1] !== null))
    ).length;
    const firstLength = files.map(f => f[0]).filter(f => f !== null).length;
    const secondLength = files.map(f => f[1]).filter(f => f !== null).length;
    if (validLength < 1) {
      return pushNotification(
        'You should select at least one input file.',
        'error'
      );
    }
    if (
      paired &&
      (firstLength !== validLength || secondLength !== validLength)
    ) {
      return pushNotification(
        'You must select the same number of mate input files.',
        'error'
      );
    }
    this.setSaving(true);
    const single = validLength === 1;
    try {
      const jobs = [];
      for (let i = 0; i < samples.length; i += 1) {
        const sample = samples[i];
        const [firstFile, secondFile] = files[i];
        if (
          firstFile !== null &&
          (!paired || (paired && secondFile !== null))
        ) {
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
                // $FlowFixMe: firstFile is not null here
                firstInputFile: firstFile.name,
                // $FlowFixMe: secondFile is not null if paired is true
                secondInputFile: paired ? secondFile.name : null
              },
              // $FlowFixMe: firstFile is not null here
              firstFile,
              paired ? secondFile : null
            )
          );
          pushNotification(`Job ${sampleName} created!`, 'success');
        }
      }
      if (!single) pushNotification('Analysis jobs created!');
      const groupJob = await this.createGroup(
        single,
        code,
        name,
        jobs,
        descriptionFile
      );
      await this.submitAnalysis(jobs, groupJob);
      refreshJobs();
      redirect(JOBS);
    } catch (e) {
      pushNotification(`An error occurred: ${e.message}`, 'error');
      this.setSaving(false);
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
              Long RNAs Analysis
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
                  length: 14
                },
                algorithm: 'salmon',
                countingAlgorithm: 'feature-counts',
                genome: 'Human_hg19_genome',
                transcriptome: 'Human_hg19_transcriptome',
                annotation: 'Human_hg19_gencode_19_gtf',
                threads: 1,
                samples: [
                  {
                    code: ''
                  }
                ]
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
                    <div>{this.getStep2(values)}</div>
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

export default withStyles(style)(LongRNA);
