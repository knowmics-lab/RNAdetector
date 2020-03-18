// @flow
import React, { Component } from 'react';
import { Typography } from '@material-ui/core';
import TextInfoContent from '@mui-treasury/components/content/textInfo';
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
        <TextInfoContent
          overline="Welcome"
          heading="RNAdetector"
          body="Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt
          ut labore et dolore magna aliqua. Orci ac auctor augue mauris."
        />
      </div>
    );
  }
}
