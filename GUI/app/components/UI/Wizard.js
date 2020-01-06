// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import FormGroup from '@material-ui/core/FormGroup';
import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import { green } from '@material-ui/core/colors';
import Stepper from '@material-ui/core/Stepper';
import Step from '@material-ui/core/Step';
import StepLabel from '@material-ui/core/StepLabel';
import { has } from 'lodash';
import { connect } from 'formik';
import type { FormikContext } from 'formik';

export type WizardProps = {
  steps: string[] | (() => string[]),
  children: React.ChildrenArray<*> | (number => React.Node),
  connectedFields?: string[][] | (number => string[]),
  formik?: FormikContext<*>,
  fieldsErrors?: { [string]: string },
  fieldsTouched?: { [string]: boolean },
  hasErrors?: number => boolean,
  prevButton?: string | (((SyntheticEvent<*>) => void) => React.Node),
  nextButton?: string | (((SyntheticEvent<*>) => void) => React.Node),
  submitButton?: string | (() => React.Node),
  onChangeActiveStep?: number => void,
  classes: {
    root: *,
    stepperContent: *,
    formControl: *,
    buttonWrapper: *,
    buttonProgress: *,
    backButton: *,
    instructions: *
  }
};

const styles = theme => ({
  root: {
    padding: theme.spacing(3, 2)
  },
  stepperContent: {
    padding: theme.spacing(3)
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
  activeStep: number
};

class Wizard extends React.Component<WizardProps, State> {
  props: WizardProps;

  static findFieldNames(component: *): string[] {
    if (!React.isValidElement(component)) return [];

    // $FlowFixMe: React.isValidElement checks that the input is a React Component
    const { props } = component;

    if (has(props, 'name')) return [props.name];
    let res = [];
    if (has(props, 'children')) {
      res = React.Children.toArray(props.children).flatMap(
        Wizard.findFieldNames
      );
    }
    return res;
  }

  static defaultProps = {
    connectedFields: [],
    fieldsErrors: [],
    fieldsTouched: [],
    formik: null,
    hasErrors: null,
    prevButton: 'Back',
    nextButton: 'Next',
    submitButton: 'Save',
    onChangeActiveStep: null
  };

  constructor(props) {
    super(props);
    this.state = {
      activeStep: 0
    };
  }

  connectedFields = null;

  handleNext = e => {
    this.setState(prev => ({
      activeStep: prev.activeStep + 1
    }));
    const { onChangeActiveStep } = this.props;
    const { activeStep } = this.state;
    if (onChangeActiveStep && typeof onChangeActiveStep === 'function') {
      onChangeActiveStep(activeStep);
    }
    e.preventDefault();
  };

  handleBack = e => {
    this.setState(prev => ({
      activeStep: prev.activeStep - 1
    }));
    const { onChangeActiveStep } = this.props;
    const { activeStep } = this.state;
    if (onChangeActiveStep && typeof onChangeActiveStep === 'function') {
      onChangeActiveStep(activeStep);
    }
    e.preventDefault();
  };

  getSteps(): string[] {
    const { steps } = this.props;
    if (typeof steps === 'function') {
      return steps();
    }
    return steps;
  }

  findConnectedFields(index: number) {
    if (this.connectedFields === null) {
      const steps = this.getSteps();
      this.connectedFields = [];
      for (let i = 0, l = steps.length; i < l; i += 1) {
        this.connectedFields.push(
          Wizard.findFieldNames(this.getStepContent(i))
        );
      }
    }
    return this.connectedFields[index] || [];
  }

  getConnectedFields(index: number): ?(string[]) {
    const { connectedFields } = this.props;
    if (!connectedFields) return null;
    if (typeof connectedFields === 'function') {
      return connectedFields(index);
    }
    return connectedFields[index] || this.findConnectedFields(index);
  }

  internalHasErrors = (connectedFields, errors, touched) => {
    if (!errors || !touched) return false;
    if (!connectedFields) return false;
    return connectedFields
      .map(
        f => has(touched, f) && touched[f] && has(errors, f) && errors[f] !== ''
      )
      .reduce((a, v) => a || v, false);
  };

  hasErrors(index: number): boolean {
    const {
      hasErrors,
      fieldsErrors: errors,
      fieldsTouched: touched,
      formik
    } = this.props;
    if (hasErrors && typeof hasErrors === 'function') {
      return hasErrors(index);
    }
    const connectedFields = this.getConnectedFields(index);
    if (formik) {
      return this.internalHasErrors(
        connectedFields,
        formik.errors,
        formik.touched
      );
    }
    return this.internalHasErrors(connectedFields, errors, touched);
  }

  getBackButton(): React.Node {
    const { prevButton, classes } = this.props;
    const { activeStep } = this.state;
    if (typeof prevButton === 'function') {
      return prevButton(this.handleBack);
    }
    const text = prevButton || 'Previous';
    return (
      <Button
        disabled={activeStep === 0}
        onClick={this.handleBack}
        className={classes.backButton}
      >
        {text}
      </Button>
    );
  }

  getNextButton(): React.Node {
    const { nextButton } = this.props;
    if (typeof nextButton === 'function') {
      return nextButton(this.handleNext);
    }
    const text = nextButton || 'Next';
    return (
      <Button variant="contained" color="primary" onClick={this.handleNext}>
        {text}
      </Button>
    );
  }

  getSubmitButton(): React.Node {
    const { submitButton } = this.props;
    if (typeof submitButton === 'function') {
      return submitButton();
    }
    const text = submitButton || 'Submit';
    return (
      <Button
        type="submit"
        variant="contained"
        color="primary"
        onClick={this.handleNext}
      >
        {text}
      </Button>
    );
  }

  getStepContent(index): React.Node {
    const { children } = this.props;
    if (typeof children === 'function') {
      return children(index);
    }
    return React.Children.toArray(children)[index] || null;
  }

  getBottomNavigator() {
    const { classes } = this.props;
    const { activeStep } = this.state;
    const steps = this.getSteps();
    return (
      <FormGroup row className={classes.formControl}>
        <Grid container justify="flex-start">
          <Grid item xs="auto">
            <div className={classes.buttonWrapper}>{this.getBackButton()}</div>
          </Grid>
          <Grid item xs="auto">
            {activeStep === steps.length - 1 ? (
              <div className={classes.buttonWrapper}>
                {this.getSubmitButton()}
              </div>
            ) : (
              <div className={classes.buttonWrapper}>
                {this.getNextButton()}
              </div>
            )}
          </Grid>
        </Grid>
      </FormGroup>
    );
  }

  render() {
    const { classes } = this.props;
    const { activeStep } = this.state;
    const steps = this.getSteps();
    return (
      <>
        <Stepper activeStep={activeStep} alternativeLabel>
          {steps.map((label, i) => (
            <Step key={label}>
              <StepLabel error={this.hasErrors(i)}>{label}</StepLabel>
            </Step>
          ))}
        </Stepper>
        <div className={classes.stepperContent}>
          {this.getStepContent(activeStep)}
          {this.getBottomNavigator()}
        </div>
      </>
    );
  }
}

Wizard.defaultProps = {
  connectedFields: [],
  fieldsErrors: [],
  fieldsTouched: [],
  formik: null,
  hasErrors: null,
  prevButton: 'Back',
  nextButton: 'Next',
  submitButton: 'Save',
  onChangeActiveStep: null
};

export const StyledWizard = withStyles(styles)(Wizard);

export default connect(StyledWizard);
