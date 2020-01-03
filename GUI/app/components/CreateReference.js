// @flow
import React, { Component } from 'react';
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

class CreateReference extends Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.uppy = Api.Uppy.initUppyInstance(['.fa', '.fasta']);
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
      availableFor: Yup.array()
        .of(Yup.string().oneOf(['bwa', 'tophat', 'hisat', 'salmon']))
        .required()
    });

  getSteps = () => ['Choose a name', 'Select aligners', 'Select a file'];

  getConnectedFields = index => {
    switch (index) {
      case 0:
        return ['name'];
      case 1:
        return ['availableFor'];
      default:
        return [];
    }
  };

  getStep0 = () => {
    const { classes } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Choose a name fo the new reference sequence. The name should contain
          only letters, numbers, and dashes.
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
          Choose which alignment algorithms will be available for this sequence.
          Based on the selected algorithms, the appropriate indexing methods
          will be used.
        </Typography>
        <SelectField
          label="Indexing algorithms"
          name="availableFor"
          options={{
            bwa: 'BWA',
            tophat: 'TopHat 2',
            hisat: 'HISAT2',
            salmon: 'Salmon'
          }}
          multiple
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
          Select the FASTA file you wish to use and click &quot;Save&quot; to
          start the upload process.
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
        const data = await Api.References.create(
          values.name,
          filename,
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
          if (await Api.Uppy.upload(this.uppy, url, pushNotification)) {
            pushNotification('FASTA file uploaded! Starting indexing job...');
            await Api.Jobs.submitJob(job.id);
            pushNotification('Indexing job queued!');
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
      pushNotification('You must select a FASTA file.', 'error');
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
            Add reference genome/transcriptome
          </Typography>
          <Formik
            initialValues={{
              name: '',
              availableFor: []
            }}
            initialErrors={validationErrors}
            validationSchema={this.getValidationSchema()}
            onSubmit={v => {
              this.formSubmit(v).catch(() => false);
            }}
          >
            {({ errors, touched }) => (
              <Form>
                <Wizard
                  fieldsErrors={errors}
                  fieldsTouched={touched}
                  connectedFields={this.getConnectedFields}
                  steps={steps}
                  submitButton={this.getSubmitButton}
                >
                  <div>{this.getStep0()}</div>
                  <div>{this.getStep1()}</div>
                  <div>{this.getStep2()}</div>
                </Wizard>
              </Form>
            )}
          </Formik>
        </Paper>
      </Box>
    );
  }
}

export default withStyles(style)(CreateReference);
