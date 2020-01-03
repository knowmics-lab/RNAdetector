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
import Stepper from '@material-ui/core/Stepper';
import Step from '@material-ui/core/Step';
import StepLabel from '@material-ui/core/StepLabel';
import Uppy from '@uppy/core';
import Tus from '@uppy/tus';
import { Dashboard } from '@uppy/react';
import { has } from 'lodash';
import * as Api from '../api';
import { JOBS } from '../constants/routes';
import SelectField from './Form/SelectField';
import TextField from './Form/TextField';

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
  activeStep: number,
  validationErrors: *
};

class CreateReference extends Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.uppy = Uppy({
      restrictions: {
        maxNumberOfFiles: 1,
        allowedFileTypes: ['.fa', '.fasta']
      },
      allowMultipleUploads: false,
      logger: Uppy.debugLogger,
      autoProceed: false
    });
    this.state = {
      isSaving: false,
      activeStep: 0,
      validationErrors: {}
    };
  }

  uppy;

  handleNext = e => {
    this.setState(prev => ({
      activeStep: prev.activeStep + 1
    }));
    e.preventDefault();
  };

  handleBack = e => {
    this.setState(prev => ({
      activeStep: prev.activeStep - 1
    }));
    e.preventDefault();
  };

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

  hasErrors = (index, errors, touched) => {
    return this.getConnectedFields(index)
      .map(
        f => has(touched, f) && touched[f] && has(errors, f) && errors[f] !== ''
      )
      .reduce((a, v) => a || v, false);
  };

  getStep0() {
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
  }

  getStep1() {
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
  }

  getStep2() {
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
                uppy={this.uppy}
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
  }

  getStepContent(index) {
    switch (index) {
      case 0:
        return this.getStep0();
      case 1:
        return this.getStep1();
      case 2:
        return this.getStep2();
      default:
        return 'Unknown step';
    }
  }

  getBottomNavigator() {
    const { classes } = this.props;
    const { isSaving, activeStep } = this.state;
    const steps = this.getSteps();
    return (
      <FormGroup row className={classes.formControl}>
        <Grid container justify="flex-start">
          <Grid item xs="auto">
            <div className={classes.buttonWrapper}>
              <Button
                disabled={activeStep === 0}
                onClick={this.handleBack}
                className={classes.backButton}
              >
                Previous
              </Button>
            </div>
          </Grid>
          <Grid item xs="auto">
            {activeStep === steps.length - 1 ? (
              <div className={classes.buttonWrapper}>
                <Button
                  type="submit"
                  variant="contained"
                  color="primary"
                  disabled={isSaving}
                >
                  Save
                </Button>
                {isSaving && (
                  <CircularProgress
                    size={24}
                    className={classes.buttonProgress}
                  />
                )}
              </div>
            ) : (
              <div className={classes.buttonWrapper}>
                <Button
                  variant="contained"
                  color="primary"
                  onClick={this.handleNext}
                >
                  Next
                </Button>
              </div>
            )}
          </Grid>
        </Grid>
      </FormGroup>
    );
  }

  initUppy = url => {
    if (this.uppy.getPlugin('Tus') !== null) {
      this.uppy.removePlugin(this.uppy.getPlugin('Tus'));
    }
    this.uppy.use(Tus, {
      endpoint: url,
      headers: {
        ...Api.Settings.getAuthHeaders()
      }
    });
  };

  checkUppyResult = (uploadResult, filename) => {
    const { pushNotification } = this.props;
    if (
      uploadResult.successful &&
      Array.isArray(uploadResult.successful) &&
      uploadResult.successful.length >= 1
    ) {
      if (uploadResult.successful[0].name === filename) {
        return true;
      }
      pushNotification(
        'Error during upload: the uploaded file is different from the selected one.',
        'error'
      );
    } else {
      pushNotification('Unknown error during upload!', 'error');
    }
    return false;
  };

  formSubmit = async values => {
    const { pushNotification, redirect, refreshJobs } = this.props;
    const file = this.uppy.getFiles()[0];
    if (file) {
      const filename = file.data.name;
      this.setState({
        isSaving: true
      });
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
          this.setState({
            isSaving: false,
            validationErrors: data.validationErrors
          });
        } else {
          const { data: job } = data;
          pushNotification(
            'A new indexing job has been created! Uploading FASTA file...'
          );
          this.initUppy(Api.Jobs.getUploadUrl(job));
          const uploadResult = await this.uppy.upload();
          if (this.checkUppyResult(uploadResult, filename)) {
            pushNotification('FASTA file uploaded! Starting indexing job...');
            await Api.Jobs.submitJob(job.id);
            pushNotification('Indexing job queued!');
            refreshJobs();
            redirect(JOBS);
          }
          this.setState({
            isSaving: false,
            validationErrors: {}
          });
        }
      } catch (e) {
        pushNotification(`An error occurred: ${e.message}`, 'error');
      }
    } else {
      pushNotification('You must select a FASTA file.', 'error');
    }
  };

  render() {
    const { classes } = this.props;
    const { activeStep, validationErrors } = this.state;
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
                <Stepper activeStep={activeStep} alternativeLabel>
                  {steps.map((label, i) => (
                    <Step key={label}>
                      <StepLabel error={this.hasErrors(i, errors, touched)}>
                        {label}
                      </StepLabel>
                    </Step>
                  ))}
                </Stepper>
                <div>
                  {this.getStepContent(activeStep)}
                  {this.getBottomNavigator()}
                </div>
              </Form>
            )}
          </Formik>
        </Paper>
      </Box>
    );
  }
}

export default withStyles(style)(CreateReference);
