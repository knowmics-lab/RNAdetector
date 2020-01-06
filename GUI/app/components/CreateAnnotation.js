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
import { Dashboard } from '@uppy/react';
import * as Api from '../api';
import { JOBS } from '../constants/routes';
import SelectField from './Form/SelectField';
import TextField from './Form/TextField';
import Wizard from './UI/Wizard';

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
  validationErrors: *
};

class CreateAnnotation extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.uppy = Api.Uppy.initUppyInstance(['.bed', '.gtf', '.gff']);
    this.state = {
      isSaving: false,
      validationErrors: {}
    };
  }

  uppy;

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

  getStep2 = () => {
    const { classes } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Select the file you wish to use and click &quot;Save&quot; to start
          the upload process.
        </Typography>
        <FormGroup row className={classes.formControl}>
          <Grid container justify="center" alignItems="center">
            <Grid item xs={12}>
              <Dashboard
                uppy={Api.Uppy.getInstance(this.uppy)}
                hideUploadButton
                proudlyDisplayPoweredByUppy={false}
                hideRetryButton
                disableStatusBar
                showLinkToFileUploadResult={false}
                width="100%"
                height="100%"
              />
            </Grid>
          </Grid>
        </FormGroup>
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
    const { pushNotification, redirect, refreshJobs } = this.props;
    if (Api.Uppy.isValid(this.uppy, pushNotification)) {
      const filename = Api.Uppy.getFilename(this.uppy, 0);
      this.setSaving(true);
      try {
        const data = await Api.Annotations.create(
          values.name,
          values.type,
          filename
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
          if (await Api.Uppy.upload(this.uppy, url, pushNotification)) {
            pushNotification('Annotation file uploaded! Starting job...');
            await Api.Jobs.submitJob(job.id);
            pushNotification('Job queued!');
            refreshJobs();
            Api.Uppy.clearInstance(this.uppy);
            this.setSaving(false);
            redirect(JOBS);
          } else {
            this.setSaving(false);
          }
        }
      } catch (e) {
        pushNotification(`An error occurred: ${e.message}`, 'error');
        this.setSaving(false);
      }
    } else {
      pushNotification('You must select a file.', 'error');
    }
  };

  render() {
    const { classes } = this.props;
    const { validationErrors } = this.state;
    const steps = this.getSteps();
    return (
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
            <Form>
              <Wizard steps={steps} submitButton={this.getSubmitButton}>
                <div>{this.getStep0()}</div>
                <div>{this.getStep1()}</div>
                <div>{this.getStep2()}</div>
              </Wizard>
            </Form>
          </Formik>
        </Paper>
      </Box>
    );
  }
}

export default withStyles(style)(CreateAnnotation);
