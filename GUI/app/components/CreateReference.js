// @flow
import React, { Component } from 'react';
import { withStyles } from '@material-ui/core/styles';
import {
  Typography,
  Paper,
  Box,
  FormGroup,
  Button,
  Grid,
  CircularProgress
} from '@material-ui/core';
import { green } from '@material-ui/core/colors';
import { Formik, Form } from 'formik';
import * as Yup from 'yup';
import Stepper from '@material-ui/core/Stepper';
import Step from '@material-ui/core/Step';
import StepLabel from '@material-ui/core/StepLabel';
import Uppy from '@uppy/core';
import Tus from '@uppy/tus';
import { DragDrop, Dashboard, StatusBar } from '@uppy/react';
import * as Api from '../api';

import SelectField from './Form/SelectField';
import TextField from './Form/TextField';

type Props = {
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
  activeStep: number
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
      activeStep: 0
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

  getSteps = () => ['Choose a name', 'Choose aligners'];

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

  getStepFinal() {
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
        return this.getStepFinal();
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
            {activeStep === steps.length ? (
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

  formSubmit = values => {
    const file = this.uppy.getFiles()[0];
    if (file) {
      const filename = file.data.name;
      Api.References.create(values.name, filename, values.availableFor)
        .then(data => {
          this.uppy.use(Tus, {
            endpoint: Api.Jobs.getUploadUrl(data),
            headers: {
              ...Api.Settings.getAuthHeaders()
            }
          });
          this.uppy
            .upload()
            .then(result => {
              console.log(result);
              return result;
            })
            .catch(e => console.log(e));
          return data;
        })
        .catch(e => console.log(e));
    }
  };

  render() {
    const { classes } = this.props;
    const { activeStep } = this.state;
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
            validationSchema={this.getValidationSchema()}
            onSubmit={this.formSubmit}
          >
            {({ values }) => (
              <Form>
                <Stepper activeStep={activeStep} alternativeLabel>
                  {steps.map(label => (
                    <Step key={label}>
                      <StepLabel>{label}</StepLabel>
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
