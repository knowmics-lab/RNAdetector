// @flow
import React, { Component } from 'react';
import { withStyles } from '@material-ui/core/styles';
import {
  Typography,
  Paper,
  Box,
  FormGroup,
  Button,
  Grid
} from '@material-ui/core';
import { green } from '@material-ui/core/colors';
import { Formik, Form } from 'formik';
import * as Yup from 'yup';
import Collapse from '@material-ui/core/Collapse';
import type { settingsStateType } from '../reducers/types';
import type { ConfigObjectType } from '../api';
import TextField from './Form/TextField';
import FileField from './Form/FileField';
import SwitchField from './Form/SwitchField';
import Snackbar from './UI/Snackbar';

type Props = {
  settings: settingsStateType,
  saveSettings: ConfigObjectType => *,
  resetSaved: () => void,
  classes: {
    root: *,
    success: *,
    icon: *,
    iconVariant: *,
    message: *,
    formControl: *
  }
};

const style = theme => ({
  root: {
    padding: theme.spacing(3, 2)
  },
  success: {
    backgroundColor: green[600]
  },
  icon: {
    fontSize: 20
  },
  iconVariant: {
    opacity: 0.9,
    marginRight: theme.spacing(1)
  },
  message: {
    display: 'flex',
    alignItems: 'center'
  },
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120
  }
});

class Settings extends Component<Props> {
  props: Props;

  handleClose = () => {
    const { resetSaved } = this.props;
    resetSaved();
  };

  formSubmit = values => {
    const { saveSettings } = this.props;
    saveSettings({
      webserviceUrl: values.webserviceUrl,
      local: values.local,
      jobsPath: values.jobsPath,
      containerName: values.containerName,
      apiKey: values.apiKey
    });
  };

  render() {
    const { classes, settings } = this.props;
    const {
      saved: isSuccessOpen,
      error: isErrorOpen,
      message: errorMessage
    } = settings.state;
    const validationSchema = Yup.object().shape({
      webserviceUrl: Yup.string().required(),
      local: Yup.boolean(),
      jobsPath: Yup.string().when('local', {
        is: true,
        then: Yup.string().required(),
        otherwise: Yup.string().notRequired()
      }),
      containerName: Yup.string().when('local', {
        is: true,
        then: Yup.string().required(),
        otherwise: Yup.string().notRequired()
      }),
      apiKey: Yup.string().required()
    });
    return (
      <Box>
        <Paper className={classes.root}>
          <Typography variant="h5" component="h3">
            Settings
          </Typography>
          <Typography component="p">Description goes here</Typography>
          <Formik
            initialValues={settings}
            validationSchema={validationSchema}
            onSubmit={this.formSubmit}
          >
            {({ values }) => (
              <Form>
                <TextField
                  label="Web Service URL"
                  name="webserviceUrl"
                  required
                />
                <SwitchField
                  label="Is docker installed locally?"
                  name="local"
                />
                <Collapse in={values.local}>
                  <FileField
                    label="Local docker storage path"
                    name="jobsPath"
                    dialogOptions={{ properties: ['openDirectory'] }}
                  />
                  <TextField
                    label="Docker container name"
                    name="containerName"
                  />
                </Collapse>
                <TextField label="API key" name="apiKey" />
                <FormGroup row className={classes.formControl}>
                  <Grid container justify="flex-end">
                    <Grid item xs="auto">
                      <Button
                        type="submit"
                        variant="contained"
                        color="primary"
                        disabled={isSuccessOpen || isErrorOpen}
                      >
                        Save
                      </Button>
                    </Grid>
                  </Grid>
                </FormGroup>
              </Form>
            )}
          </Formik>
        </Paper>
        <Snackbar
          message="Settings saved!"
          isOpen={isSuccessOpen}
          setClosed={this.handleClose}
          variant="success"
        />
        <Snackbar
          message={`An error occurred: ${errorMessage}!`}
          isOpen={isErrorOpen}
          setClosed={this.handleClose}
          variant="success"
        />
      </Box>
    );
  }
}

export default withStyles(style)(Settings);
