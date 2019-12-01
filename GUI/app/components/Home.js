// @flow
import React, { Component } from 'react';
import { Link } from 'react-router-dom';
import { Fabric } from 'office-ui-fabric-react/lib/Fabric';
import { Button, ButtonType } from 'office-ui-fabric-react/lib/Button';
import {
  Dialog,
  DialogType,
  DialogFooter
} from 'office-ui-fabric-react/lib/Dialog';
import routes from '../constants/routes';

type Props = {};
type State = {
  isOpen: boolean
};

export default class Home extends Component<Props, State> {
  props: Props;

  constructor(props: Props) {
    super(props);
    this.state = {
      isOpen: false
    };
  }

  open = () => this.setState({ isOpen: true });

  close = () => this.setState({ isOpen: false });

  render() {
    const { isOpen } = this.state;
    return (
      <Fabric className="App">
        <div style={{ margin: '5em' }}>
          <Button onClick={this.open}>I am a button.</Button>
        </div>
        <div>
          <Button>I am a button that does nothing!</Button>
        </div>
        <Dialog
          isOpen={isOpen}
          type={DialogType.close}
          onDismiss={this.close.bind(this)}
          title="Dialog title"
          subText="Dialog subText"
          isBlocking={false}
          closeButtonAriaLabel="Close"
        >
          <h1>Hello, World!</h1>
          <DialogFooter>
            <Button buttonType={ButtonType.primary} onClick={this.close}>
              OK
            </Button>
          </DialogFooter>
        </Dialog>
      </Fabric>
    );
  }
}
