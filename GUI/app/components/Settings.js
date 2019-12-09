// @flow
import React, { Component } from 'react';
import { Link } from 'react-router-dom';
import { Typography } from '@material-ui/core';
import { Formik } from 'formik';

type Props = {};

export default class Settings extends Component<Props> {
  props: Props;

  render() {
    return (
      <div>
        <Typography paragraph>Settings</Typography>
      </div>
    );
  }
}
