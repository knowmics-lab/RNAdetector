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
import Breadcrumbs from '@material-ui/core/Breadcrumbs';
import Link from '@material-ui/core/Link';
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
import type { AnalysisFileTypes, TrimGaloreConfig } from '../../types/analysis';

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
    instructions: *
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
  }
});

type State = {
  isSaving: boolean,
  firstFiles: File[],
  secondFiles: File[],
  validationErrors: *
} & UsesUpload;

class LongRNA extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.state = {
      isSaving: false,
      ...Api.Upload.ui.initUploadState(),
      firstFiles: [],
      secondFiles: [],
      validationErrors: {}
    };
  }

  getValidationSchema = () =>
    Yup.object().shape({
      // @TODO
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
                <Grid container justify="center" alignItems="center" spacing={1}>
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
          options={{
            salmon: 'Salmon',
            tophat: 'Tophat',
            hisat2: 'Hisat 2'
          }}
        />
        {(algorithm === 'tophat' || algorithm === 'hisat2') && (
          <SelectField
            label="Counting Algorithm"
            name="countingAlgorithm"
            options={{
              htseq: 'HT-seq',
              'feature-counts': 'Feature Counts',
              salmon: 'Salmon'
            }}
          />
        )}
      </>
    );
  };

  handleFileAdd = f =>
    this.setState(oldState => ({
      files: [...oldState.files, ...f]
    }));

  handleFileRemove = f =>
    this.setState(oldState => ({
      files: oldState.files.filter(o => o.path !== f.path)
    }));

  getStep2 = () => {
    const { classes } = this.props;
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
          <Grid container justify="center" alignItems="center">
            <FileSelector
              onFileRemove={this.handleFileRemove}
              onFileAdd={this.handleFileAdd}
              filters={[{ name: 'FASTA files', extensions: ['fa', 'fasta'] }]}
              disabled={isUploading}
            />
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

  formSubmit = async values => {
    const {
      pushNotification,
      redirect,
      refreshJobs,
      refreshReferences
    } = this.props;
    const { files } = this.state;
    if (files.length !== 1) {
      pushNotification('You must select a FASTA file.', 'error');
    } else {
      const file = files[0];
      this.setSaving(true);
      try {
        const data = await Api.References.create(
          values.name,
          file.name,
          values.availableFor
        );
        if (data.validationErrors) {
          pushNotification(
            'Errors occurred during validation of input parameters. Please review the form!',
            'warning'
          );
          this.setSaving(false, data.validationErrors);
        } else {
          const { data: job } = data;
          pushNotification(
            'A new indexing job has been created! Uploading FASTA file...'
          );
          const url = Api.Jobs.getUploadUrl(job);
          this.setState({
            isUploading: true,
            uploadFile: file.name
          });
          await Api.Upload.upload(
            url,
            file.path,
            file.name,
            file.type,
            (uploadedPercent, uploadedBytes, uploadTotal) =>
              this.setState({
                uploadedPercent,
                uploadedBytes,
                uploadTotal
              })
          );
          this.setState({
            isUploading: false
          });
          pushNotification('FASTA file uploaded! Starting indexing job...');
          await Api.Jobs.submitJob(job.id);
          pushNotification('Indexing job queued!');
          refreshJobs();
          refreshReferences();
          redirect(JOBS);
        }
      } catch (e) {
        pushNotification(`An error occurred: ${e.message}`, 'error');
        this.setSaving(false);
      }
    }
  };

  redirectReference = event => {
    const { redirect } = this.props;
    event.preventDefault();
    redirect(REFERENCES);
  };

  render() {
    const { classes } = this.props;
    const { validationErrors } = this.state;
    const steps = this.getSteps();
    return (
      <>
        <Box>
          <Paper className={classes.root}>
            <Typography variant="h5" component="h3">
              New Long RNAs Analysis
            </Typography>
            <Formik
              initialValues={{
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
                genome: '',
                transcriptome: '',
                annotation: '',
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
                    <div>{this.getStep2()}</div>
                  </Wizard>
                </Form>
              )}
            </Formik>
          </Paper>
        </Box>
      </>
    );
  }
}

export default withStyles(style)(LongRNA);
