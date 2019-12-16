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
  Collapse
} from '@material-ui/core';
import { green } from '@material-ui/core/colors';
import { Formik, Form } from 'formik';
import * as Yup from 'yup';
import type { settingsStateType } from '../reducers/types';
import type { ConfigObjectType } from '../api';
import TextField from './Form/TextField';
import SelectField from './Form/SelectField';
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
      local: values.local,
      apiProtocol: values.apiProtocol,
      apiHostname: values.apiHostname,
      apiPort: values.apiPort,
      apiPath: values.apiPath,
      dataPath: values.dataPath,
      dockerExecutablePath: values.dockerExecutablePath,
      containerName: values.containerName,
      apiKey: values.apiKey
    });
  };

  render() {
    const { classes, settings } = this.props;
    const {
      saving: isSaving,
      saved: isSuccessOpen,
      error: isErrorOpen,
      message: errorMessage
    } = settings.state;
    const validationSchema = Yup.object().shape({
      apiProtocol: Yup.string().required(),
      apiHostname: Yup.string().required(),
      apiPort: Yup.number()
        .positive()
        .integer()
        .min(1)
        .max(65535),
      apiPath: Yup.string().required(),
      local: Yup.boolean(),
      dataPath: Yup.string().when('local', {
        is: true,
        then: Yup.string().required(),
        otherwise: Yup.string().notRequired()
      }),
      containerName: Yup.string().when('local', {
        is: true,
        then: Yup.string().required(),
        otherwise: Yup.string().notRequired()
      }),
      apiKey: Yup.string().when('local', {
        is: true,
        then: Yup.string().notRequired(),
        otherwise: Yup.string().required()
      })
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
                <SelectField
                  label="API Protocol"
                  name="apiProtocol"
                  options={{ http: 'http', https: 'https' }}
                  required
                />
                <TextField label="API Hostname" name="apiHostname" required />
                <TextField
                  label="API Port"
                  name="apiPort"
                  type="number"
                  required
                />
                <TextField label="API Path" name="apiPath" required />
                <SwitchField
                  label="Is docker installed locally?"
                  name="local"
                />
                <Collapse in={values.local}>
                  <FileField
                    label="Local container storage path"
                    name="dataPath"
                    dialogOptions={{ properties: ['openDirectory'] }}
                  />
                  <TextField
                    label="Local container name"
                    name="containerName"
                  />
                  <FileField
                    label="Local docker executable"
                    name="dockerExecutablePath"
                    dialogOptions={{ properties: ['openFile'], filters: [] }}
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
                        disabled={isSaving}
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
          variant="error"
        />
      </Box>
    );
  }
}

export default withStyles(style)(Settings);
