import * as React from 'react';
import { ipcRenderer } from 'electron';
import Backdrop from '@material-ui/core/Backdrop';
import Typography from '@material-ui/core/Typography';
import CircularProgress from '@material-ui/core/CircularProgress';
import { makeStyles } from '@material-ui/core/styles';
import Paper from '@material-ui/core/Paper';
import Grid from '@material-ui/core/Grid';
import ErrorIcon from '@material-ui/icons/Error';

const useStyles = makeStyles(theme => ({
  backdrop: {
    zIndex: theme.zIndex.drawer + 10,
    color: '#fff'
  },
  paper: {
    padding: 10
  }
}));

export default function CloseHandler() {
  const classes = useStyles();
  const [state, setState] = React.useState({
    message: '',
    error: false,
    waiting: false
  });

  ipcRenderer.on('on-display-blocking-message', (evt, { message, error }) =>
    setState({
      message,
      error,
      waiting: true
    })
  );

  ipcRenderer.on('on-hide-blocking-message', () =>
    setState({
      message: '',
      error: false,
      waiting: false
    })
  );

  return (
    <>
      <Backdrop className={classes.backdrop} open={state.waiting}>
        <Paper elevation={3} className={classes.paper}>
          <Grid container direction="column" alignItems="center">
            {state.error ? (
              <ErrorIcon fontSize="large" color="secondary" />
            ) : (
              <CircularProgress color="inherit" />
            )}
            <Typography
              variant="h6"
              component="div"
              color={state.error ? 'secondary' : 'inherit'}
            >
              {state.message}
            </Typography>
          </Grid>
        </Paper>
      </Backdrop>
    </>
  );
}
