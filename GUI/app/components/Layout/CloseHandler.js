import * as React from 'react';
import { ipcRenderer } from 'electron';
import Backdrop from '@material-ui/core/Backdrop';
import Typography from '@material-ui/core/Typography';
import CircularProgress from '@material-ui/core/CircularProgress';
import { makeStyles } from '@material-ui/core/styles';

const useStyles = makeStyles(theme => ({
  backdrop: {
    zIndex: theme.zIndex.drawer + 10,
    color: '#fff'
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
        <CircularProgress color="inherit" />
        &nbsp;&nbsp;
        <Typography variant="h6" component="div">
          Waiting for docker to stop
        </Typography>
      </Backdrop>
    </>
  );
}
