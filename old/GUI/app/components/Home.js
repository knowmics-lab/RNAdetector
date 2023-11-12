// @flow
import React, { Component } from 'react';
import TextInfoContent from '@mui-treasury/components/content/textInfo';
import styled from 'styled-components';
import { ipcRenderer } from 'electron';
import Grid from '@material-ui/core/Grid';
import UNICT_LOGO from '../resources/unict.png';
import NERVIANO_LOGO from '../resources/nerviano.png';
import OHIO_LOGO from '../resources/ohio-logo.png';
import IOR_LOGO from '../resources/ior-logo.png';
import { GUI_VERSION } from '../constants/system.json';
import { Settings, Utils } from '../api';
import type { Capabilities } from '../types/common';

type Props = {};

type State = {
  capabilities: ?Capabilities
};

const FooterContainer = styled.div`
  text-align: center;
  margin-top: 200px;
`;

export default class Home extends Component<Props, State> {
  props: Props;

  constructor(props: Props) {
    super(props);
    this.state = {
      capabilities: Utils.capabilities()
    };
  }

  componentDidMount() {
    if (!Utils.capabilitiesLoaded()) {
      if (!Settings.isLocal()) {
        ipcRenderer.send('display-blocking-message', {
          message: 'Loading...',
          error: false
        });
      }
      let timer;
      const doRefresh = async () => {
        if ((await this.refreshCapabilities()) && timer) {
          if (!Settings.isLocal()) {
            ipcRenderer.send('hide-blocking-message');
          }
          clearInterval(timer);
        }
      };
      timer = setInterval(doRefresh, 2000);
    }
  }

  async refreshCapabilities() {
    try {
      const capabilities = await Utils.refreshCapabilities();
      this.setState({
        capabilities
      });
      return true;
    } catch (e) {
      return false;
    }
  }

  render() {
    const { capabilities } = this.state;
    return (
      <div>
        <TextInfoContent
          overline="Welcome"
          heading="RNAdetector"
          body="Welcome to RNAdetector the user-friendly open-source pipeline for RNA-seq data analysis."
        />
        <TextInfoContent
          body={`You are currently using version ${GUI_VERSION} which uses the RNAdetector Docker
          Container v. ${
            capabilities
              ? capabilities.containerVersion
              : '(Loading...Please Wait...)'
          }`}
        />
        <FooterContainer>
          <Grid container alignItems="center" justify="center" spacing={2}>
            <Grid item xs={12} sm={6} md={4} xl={3}>
              <img
                src={UNICT_LOGO}
                alt="UNICT"
                style={{ height: '50px', width: 'auto' }}
              />
            </Grid>
            <Grid item xs={12} sm={6} md={4} xl={3}>
              <img
                src={NERVIANO_LOGO}
                alt="NMS"
                style={{ height: '50px', width: 'auto' }}
              />
            </Grid>
            <Grid item xs={12} sm={6} md={4} xl={3}>
              <img
                src={IOR_LOGO}
                alt="IOR"
                style={{ height: '50px', width: 'auto' }}
              />
            </Grid>
            <Grid item xs={12} sm={6} md={4} xl={3}>
              <img
                src={OHIO_LOGO}
                alt="OSU"
                style={{ height: '50px', width: 'auto' }}
              />
            </Grid>
          </Grid>
        </FooterContainer>
      </div>
    );
  }
}
