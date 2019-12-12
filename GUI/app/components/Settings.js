// @flow
import React, { Component } from 'react';
import { withStyles } from '@material-ui/core/styles';
import { Link } from 'react-router-dom';
import { Typography, Paper, Box } from '@material-ui/core';
import { Formik, Form, Field } from 'formik';
import * as Yup from 'yup';
import type { settingsStateType } from '../reducers/types';
import TextField from './Form/TextField';
import SelectField from './Form/SelectField';
import FileField from './Form/FileField';

type Props = {
  settings: settingsStateType,
  classes: {
    root: *
  }
};

const style = theme => ({
  root: {
    padding: theme.spacing(3, 2)
  }
});

class Settings extends Component<Props> {
  props: Props;

  formSubmit = (values, helpers) => {
    console.log(values);
  };

  render() {
    const { classes, settings } = this.props;
    const validationSchema = Yup.object().shape({});
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
            {({ errors, touched, isValidating }) => (
              <Form>
                <TextField
                  label="Web Service URL"
                  name="webserviceUrl"
                  required
                />
                <FileField
                  label="Web Service URL"
                  name="webserviceUrl"
                  required
                />
              </Form>
            )}
          </Formik>
        </Paper>
      </Box>
    );
  }
}

export default withStyles(style)(Settings);
