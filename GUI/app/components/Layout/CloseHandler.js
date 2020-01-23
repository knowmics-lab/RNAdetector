import * as React from 'react';
import { ipcRenderer } from 'electron';
import Backdrop from '@material-ui/core/Backdrop';
import Typography from '@material-ui/core/Typography';
import CircularProgress from '@material-ui/core/CircularProgress';
import { makeStyles } from '@material-ui/core/styles';
import Paper from '@material-ui/core/Paper';
import Grid from '@material-ui/core/Grid';

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
  const [waiting, setWaiting] = React.useState(false);

  ipcRenderer.on('waiting-docker-close', () => {
    setWaiting(true);
  });

  return (
    <>
      <Backdrop className={classes.backdrop} open={waiting}>
        <Paper elevation={3} className={classes.paper}>
          <Grid container direction="column" alignItems="center">
            <CircularProgress color="inherit" />
            &nbsp;&nbsp;
            <Typography variant="h6" component="div">
              Waiting for docker to stop
            </Typography>
          </Grid>
        </Paper>
      </Backdrop>
    </>
  );
}
