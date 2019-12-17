// @flow
import React, { Component } from 'react';
import { Typography } from '@material-ui/core';
import Table from './UI/Table';

type Props = {};

export default class Home extends Component<Props> {
  props: Props;

  render() {
    return (
      <div>
        <Typography paragraph>Hello World</Typography>
        <Table columns={[

        ]}/>
      </div>
    );
  }
}
