// @flow
import React, { Component } from 'react';
import TextInfoContent from '@mui-treasury/components/content/textInfo';

type Props = {};

export default class Home extends Component<Props> {
  props: Props;

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
