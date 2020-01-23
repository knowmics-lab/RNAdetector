import * as React from 'react';
import Backdrop from '@material-ui/core/Backdrop';
import Typography from '@material-ui/core/Typography';
import CircularProgress from '@material-ui/core/CircularProgress';
import { makeStyles } from '@material-ui/core/styles';
import Paper from '@material-ui/core/Paper';
import Grid from '@material-ui/core/Grid';
import ErrorIcon from '@material-ui/icons/Error';
import * as Api from '../../api';

const useStyles = makeStyles(theme => ({
  backdrop: {
    zIndex: theme.zIndex.drawer + 10,
    color: '#fff'
  },
  paper: {
    padding: 10
  }
}));

export default function StartHandler() {
  const classes = useStyles();
  const [first, setFirst] = React.useState(false);
  const [state, setState] = React.useState({
    waiting: false,
    message: ''
  });

  React.useEffect(() => {
    if (!first) {
      if (Api.Settings.isConfigured() && Api.Settings.isLocal()) {
        setState({
          waiting: true,
          error: false,
          message: 'Starting docker container...'
        });
        Api.Docker.startContainer()
          // eslint-disable-next-line promise/always-return
          .then(() => {
            setState({
              waiting: false,
              error: false,
              message: ''
            });
          })
          .catch(e => {
            setState({
              waiting: true,
              error: true,
              message: e.message
            });
          });
      }
      setFirst(true);
    }
  });

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
