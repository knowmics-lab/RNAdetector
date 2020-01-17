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
import * as Api from '../../api';
import { JOBS, REFERENCES } from '../../constants/routes';
import SelectField from '../Form/SelectField';
import TextField from '../Form/TextField';
import Wizard from '../UI/Wizard';
import FileSelector from '../UI/FileSelector';
import type { File } from '../UI/FileSelector';
import UploadProgress from '../UI/UploadProgress';
import type { UsesUpload } from '../../types/ui';
import SwitchField from '../Form/SwitchField';
import type {
  AnalysisFileTypes,
  LongRNAAnalysisConfig,
  TrimGaloreConfig
} from '../../types/analysis';
import type { SimpleMapType } from '../../types/common';

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
  isSaving: boolean,
  firstFiles: File[],
  secondFiles: File[],
  genomes: SimpleMapType<SimpleMapType<string>>,
  annotations: SimpleMapType<string>,
  validationErrors: *
} & UsesUpload;

class LongRNA extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.state = {
      isLoading: false,
      isSaving: false,
      ...Api.Upload.ui.initUploadState(),
      firstFiles: [],
      secondFiles: [],
      validationErrors: {},
      genomes: {},
      annotations: {}
    };
  }

  componentDidMount(): void {
    const { pushNotification } = this.props;
    Object.keys(ALGORITHMS).forEach((a: string) => {
      this.setState({ isLoading: true });
      Api.References.fetchAllByAlgorithm(a === 'hisat2' ? 'hisat' : a)
        .then(v =>
          this.setState(prevState => ({
            isLoading: false,
            genomes: { ...prevState.genomes, [a]: v }
          }))
        )
        .catch(e => {
          pushNotification(`An error occurred: ${e.message}!`, 'error');
          this.setState({ isLoading: false });
        });
    });
    this.setState({ isLoading: true });
    Api.Annotations.fetchAllByType('gtf')
      .then(v =>
        this.setState({
          isLoading: false,
          annotations: v
        })
      )
      .catch(e => {
        pushNotification(`An error occurred: ${e.message}!`, 'error');
        this.setState({ isLoading: false });
      });
  }

  getValidationSchema = () =>
    Yup.object().shape({
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
          Choose the type of input files and sequencing strategy (single or
          paired-end).
        </Typography>
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
          Choose which steps will be included in the analysis: trimming, BAM to
          FASTQ conversion, alignment and counting, or quantification.
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
          />
        )}
        {(algorithm === 'salmon' || countingAlgorithm === 'salmon') && (
          <SelectField
            label="Reference Transcriptome"
            name="transcriptome"
            options={genomes.salmon}
          />
        )}
        {((algorithm !== 'salmon' && countingAlgorithm !== 'salmon') ||
          algorithm === 'tophat') && (
          <SelectField
            label="Genome Annotation"
            name="annotation"
            options={annotations}
          />
        )}
      </>
    );
  };

  makeHandleFileAdd = (field: string) => f =>
    this.setState(oldState => ({
      [field]: [...oldState[field], ...f]
    }));

  makeHandleFileRemove = (field: string) => f =>
    this.setState(oldState => ({
      [field]: oldState[field].filter(o => o.path !== f.path)
    }));

  getStep3 = values => {
    const { classes } = this.props;
    const { inputType, paired } = values;
    const {
      isUploading,
      uploadFile,
      uploadedBytes,
      uploadedPercent,
      uploadTotal
    } = this.state;
    return (
      <>
        <Typography className={classes.instructions}>
          Select the FASTA file you wish to use and click &quot;Save&quot; to
          start the upload process.
        </Typography>
        <FormGroup row className={classes.formControl}>
          <Grid container justify="space-around" alignItems="flex-start">
            <FileSelector
              onFileRemove={this.makeHandleFileRemove('firstFiles')}
              onFileAdd={this.makeHandleFileAdd('firstFiles')}
              filters={Api.Utils.analysisFileExtensions(inputType)}
              disabled={isUploading}
              multiple
            />
            {paired && (
              <FileSelector
                onFileRemove={this.makeHandleFileRemove('secondFiles')}
                onFileAdd={this.makeHandleFileAdd('secondFiles')}
                filters={Api.Utils.analysisFileExtensions(inputType)}
                disabled={isUploading}
                multiple
              />
            )}
          </Grid>
        </FormGroup>
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
          Save
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

  uploadFile = async (url: string, file: File) => {
    Api.Upload.ui.uploadStart(this.setState.bind(this), file.name);
    await Api.Upload.upload(
      url,
      file.path,
      file.name,
      file.type,
      Api.Upload.ui.makeOnProgress(this.setState.bind(this))
    );
    Api.Upload.ui.uploadEnd(this.setState.bind(this));
  };

  analysisSubmit = async (
    name: string,
    parameters: LongRNAAnalysisConfig,
    firstFile: File,
    secondFile: ?File
  ) => {
    const { pushNotification } = this.props;
    const data = await Api.Analysis.createLongRNA(name, parameters);
    if (data.validationErrors) {
      pushNotification(
        'Errors occurred during validation of input parameters. Please review the form!',
        'warning'
      );
      this.setSaving(false, data.validationErrors);
    } else {
      const { data: job } = data;
      pushNotification('Job created! Uploading files...');
      const url = Api.Jobs.getUploadUrl(job);
      await this.uploadFile(url, firstFile);
      if (secondFile) await this.uploadFile(url, secondFile);
      pushNotification('Files uploaded! Starting job...');
      await Api.Jobs.submitJob(job.id);
      pushNotification('Job queued!');
    }
  };

  formSubmit = async values => {
    const { paired } = values;
    const { name, ...params } = values;
    const { pushNotification, redirect, refreshJobs } = this.props;
    const { firstFiles, secondFiles } = this.state;
    if (firstFiles.length < 1) {
      return pushNotification(
        'You should select at least one input file.',
        'error'
      );
    }
    if (paired && secondFiles.length !== firstFiles.length) {
      return pushNotification(
        'You must select the same number of mate input files.',
        'error'
      );
    }
    this.setSaving(true);
    const single = firstFiles.length === 1;
    const promises = firstFiles.map((file1, i) => {
      return this.analysisSubmit(
        `${name}${single ? '' : ` - Sample ${i}`}`,
        {
          ...params,
          firstInputFile: file1.name,
          secondInputFile: paired ? secondFiles[i].name : null
        },
        file1,
        paired ? secondFiles[i] : null
      );
    });
    try {
      await Promise.all(promises);
      refreshJobs();
      redirect(JOBS);
    } catch (e) {
      console.log(e);
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
                threads: 1
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