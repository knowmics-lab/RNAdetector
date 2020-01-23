// @flow
import React, { Component } from 'react';
import { Typography } from '@material-ui/core';
import { SETTINGS } from '../constants/routes';
import * as Api from '../api';

type Props = {
  redirect: mixed => void
};

export default class Home extends Component<Props> {
  props: Props;

  componentDidMount() {
    if (!Api.Settings.isConfigured()) {
      const { redirect } = this.props;
      redirect(SETTINGS);
    }
  }

  render() {
    return (
      <div>
        <Typography paragraph>TODO</Typography>
      </div>
    );
  }
}
