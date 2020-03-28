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
import * as Api from '../api';
import { JOBS, ANNOTATIONS } from '../constants/routes';
import SelectField from './Form/SelectField';
import TextField from './Form/TextField';
import Wizard from './UI/Wizard';
import type { File } from './UI/FileSelector';
import FileSelector from './UI/FileSelector';
import UploadProgress from './UI/UploadProgress';

type Props = {
  refreshJobs: () => void,
  refreshAnnotations: () => void,
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
  isUploading: boolean,
  uploadFile: string,
  uploadedBytes: number,
  uploadedPercent: number,
  uploadTotal: number,
  annotationFile: ?File,
  mapFile: ?File,
  validationErrors: *
};

class CreateAnnotation extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.state = {
      isSaving: false,
      isUploading: false,
      uploadFile: '',
      uploadedBytes: 0,
      uploadedPercent: 0,
      uploadTotal: 0,
      annotationFile: null,
      mapFile: null,
      validationErrors: {}
    };
  }

  getValidationSchema = () =>
    Yup.object().shape({
      name: Yup.string()
        .required()
        .matches(/^[A-Za-z0-9\-_]+$/, {
          message: 'The field must contain only letters, numbers, and dashes.'
        }),
      type: Yup.string().oneOf(['gtf', 'bed'])
    });

  getSteps = () => ['Choose a name', 'Select a type', 'Select a file'];

  getStep0 = () => {
    const { classes } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Choose a name fo the new annotation. The name should contain only
          letters, numbers, and dashes.
        </Typography>
        <TextField label="Name" name="name" required />
      </>
    );
  };

  getStep1 = () => {
    const { classes } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Choose the type of annotation (GTF or BED).
        </Typography>
        <SelectField
          label="Annotation type"
          name="type"
          options={{
            gtf: 'GTF',
            bed: 'BED'
          }}
          required
        />
      </>
    );
  };

  handleFileAdd = (field: string) => f =>
    this.setState({
      [field]: f[0]
    });

  handleFileRemove = (field: string) => () =>
    this.setState({
      [field]: null
    });

  getStep2 = values => {
    const { classes } = this.props;
    const {
      isUploading,
      uploadFile,
      uploadedBytes,
      uploadedPercent,
      uploadTotal
    } = this.state;
    const filters = [];
    if (values.type === 'bed') {
      filters.push({ name: 'BED files', extensions: ['bed'] });
    } else {
      filters.push({ name: 'GTF files', extensions: ['gtf'] });
    }
    return (
      <>
        <Typography className={classes.instructions}>
          Select the annotation file and an optional map to Entrez id if you
          wish to use pathway analysis and click &quot;Save&quot; to start the
          upload process.
        </Typography>
        <FormGroup row className={classes.formControl}>
          <Grid container justify="center" alignItems="flex-start">
            <Grid item xs>
              <FileSelector
                title="Select annotation file"
                onFileRemove={this.handleFileRemove('annotationFile')}
                onFileAdd={this.handleFileAdd('annotationFile')}
                filters={filters}
                disabled={isUploading}
              />
            </Grid>
            <Grid item xs>
              <FileSelector
                title="Select map file"
                onFileRemove={this.handleFileRemove('mapFile')}
                onFileAdd={this.handleFileAdd('mapFile')}
                filters={[{ name: 'TSV files', extensions: ['tsv', 'txt'] }]}
                disabled={isUploading}
              />
            </Grid>
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

  redirectAnnotation = event => {
    const { redirect } = this.props;
    event.preventDefault();
    redirect(ANNOTATIONS);
  };

  formSubmit = async values => {
    const {
      pushNotification,
      redirect,
      refreshJobs,
      refreshAnnotations
    } = this.props;
    const { annotationFile, mapFile } = this.state;
    if (!annotationFile) {
      pushNotification('You must select an annotation file.', 'error');
    } else {
      this.setSaving(true);
      try {
        const data = await Api.Annotations.create(
          values.name,
          values.type,
          annotationFile.name,
          mapFile ? mapFile.name : null
        );
        if (data.validationErrors) {
          pushNotification(
            'Errors occurred during validation of input parameters. Please review the form!',
            'warning'
          );
          this.setSaving(false, data.validationErrors);
        } else {
          const { data: job } = data;
          pushNotification('Job created! Uploading files...');
          const files = [annotationFile];
          if (mapFile) files.push(mapFile);
          // eslint-disable-next-line no-restricted-syntax
          for (const f of files) {
            this.setState({
              isUploading: true,
              uploadFile: f.name,
              uploadedBytes: 0,
              uploadedPercent: 0,
              uploadTotal: 0
            });
            // eslint-disable-next-line no-await-in-loop
            await Api.Upload.upload(
              job,
              f.path,
              f.name,
              f.type,
              (uploadedPercent, uploadedBytes, uploadTotal) =>
                this.setState({
                  uploadedPercent,
                  uploadedBytes,
                  uploadTotal
                })
            );
          }
          this.setState({
            isUploading: false
          });
          pushNotification('Annotation file uploaded! Starting job...');
          await Api.Jobs.submitJob(job.id);
          pushNotification('Job queued!');
          refreshJobs();
          refreshAnnotations();
          redirect(JOBS);
        }
      } catch (e) {
        pushNotification(`An error occurred: ${e.message}`, 'error');
        this.setSaving(false);
      }
    }
  };

  render() {
    const { classes } = this.props;
    const { validationErrors } = this.state;
    const steps = this.getSteps();
    return (
      <>
        <Breadcrumbs aria-label="breadcrumb">
          <Link color="inherit" href="#" onClick={this.redirectAnnotation}>
            Annotations
          </Link>
          <Typography color="textPrimary">Add Annotation</Typography>
        </Breadcrumbs>
        <Box>
          <Paper className={classes.root}>
            <Typography variant="h5" component="h3">
              Add annotation
            </Typography>
            <Formik
              initialValues={{
                name: '',
                type: 'gtf'
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
                    <div>{this.getStep1()}</div>
                    <div>{this.getStep2(values)}</div>
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

export default withStyles(style)(CreateAnnotation);
