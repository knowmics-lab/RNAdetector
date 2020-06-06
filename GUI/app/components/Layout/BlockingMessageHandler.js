// eslint-disable-next-line max-classes-per-file
import * as React from 'react';
import { ipcRenderer } from 'electron';
import Backdrop from '@material-ui/core/Backdrop';
import Typography from '@material-ui/core/Typography';
import CircularProgress from '@material-ui/core/CircularProgress';
import { withStyles } from '@material-ui/core/styles';
import Paper from '@material-ui/core/Paper';
import Grid from '@material-ui/core/Grid';
import ErrorIcon from '@material-ui/icons/Error';

const style = theme => ({
  backdrop: {
    zIndex: theme.zIndex.drawer + 10,
    color: '#fff'
  },
  paper: {
    padding: 10
  }
});

type LogState = {
  log: string
};

class LogHandler extends React.Component<{}, LogState> {
  constructor(props) {
    super(props);
    this.state = {
      log: ''
    };
  }

  componentDidMount() {
    ipcRenderer.on('on-blocking-message-log', (evt, log) =>
      this.setState({
        log
      })
    );
  }

  shouldComponentUpdate(nextProps, nextState) {
    const { log: oldLog } = this.state;
    return oldLog !== nextState.log;
  }

  componentWillUnmount() {
    ipcRenderer.removeAllListeners('on-blocking-message-log');
  }

  render() {
    const { log } = this.state;
    return <>{log ? <pre>{log}</pre> : null}</>;
  }
}

type CloseProps = {
  classes: {
    backdrop: *,
    paper: *
  }
};

type CloseState = {
  message: string,
  error: boolean,
  waiting: boolean
};

class BlockingMessageHandler extends React.Component<CloseProps, CloseState> {
  props: CloseProps;

  constructor(props) {
    super(props);
    this.state = {
      message: '',
      error: false,
      waiting: false
    };
  }

  componentDidMount() {
    ipcRenderer.on('on-display-blocking-message', (evt, { message, error }) =>
      this.setState({
        message,
        error,
        waiting: true
      })
    );

    ipcRenderer.on('on-hide-blocking-message', () =>
      this.setState({
        message: '',
        error: false,
        waiting: false
      })
    );
  }

  componentWillUnmount() {
    ipcRenderer.removeAllListeners('on-display-blocking-message');
    ipcRenderer.removeAllListeners('on-hide-blocking-message');
  }

  render() {
    const { classes } = this.props;
    const { message, error, waiting } = this.state;
    return (
      <>
        <Backdrop className={classes.backdrop} open={waiting}>
          <Paper elevation={3} className={classes.paper}>
            <Grid container direction="column" alignItems="center">
              {error ? (
                <ErrorIcon fontSize="large" color="secondary" />
              ) : (
                <CircularProgress color="inherit" />
              )}
              <Typography
                variant="h6"
                component="div"
                color={error ? 'secondary' : 'inherit'}
              >
                {message}
              </Typography>
              <LogHandler />
            </Grid>
          </Paper>
        </Backdrop>
      </>
    );
  }
}

export default withStyles(style)(BlockingMessageHandler);
