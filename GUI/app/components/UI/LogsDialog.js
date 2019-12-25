/* eslint-disable no-nested-ternary */
import React from 'react';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import useMediaQuery from '@material-ui/core/useMediaQuery';
import { useTheme } from '@material-ui/core/styles';
import { TextField } from '@material-ui/core';
import InputAdornment from '@material-ui/core/InputAdornment';
import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import { has } from 'lodash';
import CircularProgress from '@material-ui/core/CircularProgress';
import Icon from '@material-ui/core/Icon';
import type { StateType } from '../../reducers/types';
import * as JobsActions from '../../actions/jobs';
import type { LoadedJobs } from '../../types/jobs';
import { ContentWrapper as Snackbar } from './Snackbar';

type LogsDialogProps = {|
  jobId: ?number,
  open: boolean,
  onClose: () => void
|};

type InternalLogsDialogProps = LogsDialogProps & {
  jobs: LoadedJobs,
  requestJob: (number, boolean) => void
};

function InternalLogsDialog({
  jobId,
  open,
  onClose,
  jobs,
  requestJob
}: InternalLogsDialogProps) {
  const theme = useTheme();
  const fnRequestJob = () => {
    requestJob(jobId, true);
  };
  const [timeout, setTimeout] = React.useState(30);
  const [timer, setTimer] = React.useState();
  const fullScreen = useMediaQuery(theme.breakpoints.down('md'));

  const isOpen = !!jobId && open;

  const hasJob = jobId && has(jobs.items, jobId);

  const startTimer = () => {
    if (timer) {
      clearInterval(timer);
    }
    setTimer(setInterval(fnRequestJob, timeout * 1000));
  };

  const stopTimer = () => {
    if (timer) {
      clearInterval(timer);
      setTimer(null);
    }
  };

  const handleTimeoutChange = e => {
    const val = +e.target.value;
    if (timeout !== val) {
      setTimeout(val);
      if (timer) {
        startTimer();
      }
    }
  };

  const needsRefresh = hasJob && jobs.items[jobId].status === 'processing';

  if (isOpen && !hasJob && !jobs.meta.fetching) {
    fnRequestJob();
  }

  if (isOpen && !timer && needsRefresh) {
    startTimer();
  }

  if (!isOpen && timer && needsRefresh) {
    stopTimer();
  }

  return (
    <>
      <Dialog fullScreen={fullScreen} open={isOpen} onClose={onClose}>
        <DialogTitle>
          {hasJob ? `Logs of ${jobs.items[jobId].name}` : 'Logs'}
        </DialogTitle>
        <DialogContent>
          {hasJob ? (
            <pre>{jobs.items[jobId].log}</pre>
          ) : jobs.meta.fetching && jobs.meta.error ? (
            <div style={{ textAlign: 'center' }}>
              <Snackbar
                message={`An error occurred: ${jobs.meta.error}!`}
                onClose={null}
                variant="error"
              />
            </div>
          ) : (
            <div style={{ textAlign: 'center' }}>
              <CircularProgress />
            </div>
          )}
        </DialogContent>
        <DialogActions>
          {timer ? (
            <>
              <Icon className="fas fa-sync fa-spin" />
              <TextField
                label="Refresh every"
                variant="filled"
                value={timeout}
                onChange={handleTimeoutChange}
                InputProps={{
                  endAdornment: (
                    <InputAdornment position="end">s</InputAdornment>
                  )
                }}
              />
            </>
          ) : (
            <Icon className="fas fa-sync" />
          )}
          <Button onClick={onClose} color="primary" autoFocus>
            Close
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
}

export default connect(
  (state: StateType) => ({
    jobs: state.jobs.jobs
  }),
  dispatch =>
    bindActionCreators(
      {
        requestJob: JobsActions.requestJob
      },
      dispatch
    )
)(InternalLogsDialog);
