// @flow
import React, { Component } from 'react';
import { Typography } from '@material-ui/core';
import { Jobs } from '../api';

type Props = {};

Jobs.fetchJobById(5).then(data => console.log(data));

export default class Home extends Component<Props> {
  props: Props;

  render() {
    return (
      <div>
        <Typography paragraph>Hello World</Typography>
      </div>
    );
  }
}
