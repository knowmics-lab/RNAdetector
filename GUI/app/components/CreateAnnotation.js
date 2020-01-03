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
import Stepper from '@material-ui/core/Stepper';
import Step from '@material-ui/core/Step';
import StepLabel from '@material-ui/core/StepLabel';
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

class CreateAnnotation extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.uppy = Api.Uppy.initUppyInstance(['.bed', '.gtf', '.gff']);
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
      type: Yup.string().oneOf(['gtf', 'bed'])
    });

  getSteps = () => ['Choose a name', 'Select a type', 'Select a file'];

  getConnectedFields = index => {
    switch (index) {
      case 0:
        return ['name'];
      case 1:
        return ['type'];
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
          Choose a name fo the new annotation. The name should contain only
          letters, numbers, and dashes.
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
  }

  getStep2() {
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

  formSubmit = async values => {
    const { pushNotification, redirect, refreshJobs } = this.props;
    if (Api.Uppy.isValid(this.uppy, pushNotification)) {
      const filename = Api.Uppy.getFilename(this.uppy, 0);
      this.setState({
        isSaving: true
      });
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
          this.setState({
            isSaving: false,
            validationErrors: data.validationErrors
          });
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
            redirect(JOBS);
            Api.Uppy.clearInstance(this.uppy);
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
      pushNotification('You must select a file.', 'error');
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
            Add annotation genome/transcriptome
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

export default withStyles(style)(CreateAnnotation);
