// @flow
import React, { Component } from 'react';
import { Link } from 'react-router-dom';
import { Typography } from '@material-ui/core';
import { Formik } from 'formik';
import * as Yup from 'yup';

type Props = {};

export default class Settings extends Component<Props> {
  props: Props;

  formSubmit = (values, helpers) => {
    console.log(values);
  }

  render() {
    const validationSchema = Yup.object().shape({

    });
    return (
      <div>
        <Formik initialValues=[{}] validationSchema={validationSchema} onSubmit={this.formSubmit}>
          {({ errors, touched, isValidating}) => {

          }}
        </Formik>
        <Typography paragraph>Settings</Typography>
      </div>
    );
  }
}
