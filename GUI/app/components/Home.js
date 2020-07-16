// @flow
import React, { Component } from 'react';
import TextInfoContent from '@mui-treasury/components/content/textInfo';
import styled from 'styled-components';
import GridList from '@material-ui/core/GridList';
import GridListTile from '@material-ui/core/GridListTile';
import UNICT_LOGO from '../resources/unict.png';
import { DOCKER_IMAGE_VERSION, GUI_VERSION } from '../constants/system.json';

type Props = {};

const FooterContainer = styled.div`
  text-align: center;
  position: fixed;
  bottom: 0;
  margin-bottom: 20px;
  width: 100% !important;
`;

export default class Home extends Component<Props> {
  props: Props;

  render() {
    return (
      <div>
        <TextInfoContent
          overline="Welcome"
          heading="RNAdetector"
          body="Welcome to RNAdetector the user-friendly open-source pipeline for RNA-seq data analysis."
        />
        <TextInfoContent
          body={`You are currently using version ${GUI_VERSION} which uses the RNAdetector Docker
          Container v. ${DOCKER_IMAGE_VERSION}`}
        />
        <FooterContainer>
          <GridList cellHeight={100} cols={3}>
            <GridListTile cols={1}>
              <img src={UNICT_LOGO} alt="UNICT" style={{ height: '100px', width: 'auto' }} />
            </GridListTile>
          </GridList>
        </FooterContainer>
      </div>
    );
  }
}
