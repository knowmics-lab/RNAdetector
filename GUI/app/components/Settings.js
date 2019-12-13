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
  Icon
} from '@material-ui/core';
import { green } from '@material-ui/core/colors';
import { Formik, Form } from 'formik';
import * as Yup from 'yup';
import Collapse from '@material-ui/core/Collapse';
import Snackbar from '@material-ui/core/Snackbar';
import SnackbarContent from '@material-ui/core/SnackbarContent';
import CloseIcon from '@material-ui/icons/Close';
import IconButton from '@material-ui/core/IconButton';
import CheckCircleIcon from '@material-ui/icons/CheckCircle';
import type { settingsStateType } from '../reducers/types';
import TextField from './Form/TextField';
import FileField from './Form/FileField';
import SwitchField from './Form/SwitchField';

type Props = {
  settings: settingsStateType,
  saveSettings: settingsStateType => *,
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

type SettingsState = {
  successOpen: boolean
};

class Settings extends Component<Props, SettingsState> {
  props: Props;

  constructor(props, context) {
    super(props, context);

    this.state = {
      successOpen: false
    };
  }

  handleOpen = () => {
    this.setState({
      successOpen: true
    });
  };

  handleClose = (event, reason) => {
    if (reason === 'clickaway') {
      return;
    }

    this.setState({
      successOpen: false
    });
  };

  formSubmit = values => {
    const { saveSettings } = this.props;
    saveSettings({
      webserviceUrl: values.webserviceUrl,
      local: values.local,
      jobsPath: values.jobsPath
    });
    this.handleOpen();
  };

  render() {
    const { successOpen } = this.state;
    const { classes, settings } = this.props;
    const validationSchema = Yup.object().shape({
      webserviceUrl: Yup.string().required(),
      local: Yup.boolean(),
      jobsPath: Yup.string().when('local', {
        is: true,
        then: Yup.string().required(),
        otherwise: Yup.string().notRequired()
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
                </Collapse>
                <FormGroup row className={classes.formControl}>
                  <Grid container justify="flex-end">
                    <Grid item xs="auto">
                      <Button type="submit" variant="contained" color="primary">
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
          anchorOrigin={{
            vertical: 'bottom',
            horizontal: 'right'
          }}
          open={successOpen}
          autoHideDuration={3000}
          onClose={this.handleClose}
        >
          <SnackbarContent
            className={classes.success}
            message={
              <span className={classes.message}>
                <CheckCircleIcon />
                Settings saved!
              </span>
            }
            action={[
              <IconButton
                key="close"
                aria-label="close"
                color="inherit"
                onClick={this.handleClose}
              >
                <CloseIcon className={classes.icon} />
              </IconButton>
            ]}
          />
        </Snackbar>
      </Box>
    );
  }
}

export default withStyles(style)(Settings);
