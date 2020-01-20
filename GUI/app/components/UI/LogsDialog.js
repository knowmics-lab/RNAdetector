// @flow
/* eslint-disable no-nested-ternary */
import * as React from 'react';
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

type InternalLogsDialogProps = {
  jobs: LoadedJobs,
  requestJob: (number, boolean) => void,
  jobId: ?number,
  open: boolean,
  onClose: () => void
};

function InternalLogsDialog({
  jobId,
  open,
  onClose,
  jobs,
  requestJob
}: InternalLogsDialogProps) {
  const logRef = React.createRef();
  const theme = useTheme();
  const fnRequestJob = () => {
    if (jobId) requestJob(jobId, true);
  };
  const [first, setFirst] = React.useState(true);
  const [timeout, setTimeout] = React.useState(30);
  const [timer, setTimer] = React.useState();
  const fullScreen = useMediaQuery(theme.breakpoints.down('md'));
  const isOpen = !!jobId && open;
  const hasJob = jobId && has(jobs.items, jobId);
  const startTimer = (val: ?number = null) => {
    if (timer) {
      clearInterval(timer);
    }
    setTimer(setInterval(fnRequestJob, (val || timeout) * 1000));
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
      if (timer) startTimer(val);
    }
  };
  const internalOnClose = () => {
    stopTimer();
    setFirst(true);
    onClose();
  };
  const needsRefresh =
    jobId && hasJob && jobs.items[jobId].status === 'processing';
  const doRequestJob = isOpen && !hasJob && !jobs.fetching;
  const doStartTimer = isOpen && !timer && needsRefresh;
  const doStopTimer =
    (!isOpen && timer && needsRefresh) || (isOpen && timer && !needsRefresh);
  if (doRequestJob) fnRequestJob();
  if (doStartTimer) startTimer();
  if (doStopTimer) stopTimer();
  if (isOpen && first) {
    fnRequestJob();
    setFirst(false);
  }

  React.useEffect(() => {
    if (logRef.current) {
      logRef.current.scrollIntoView({ behavior: 'smooth' });
    }
  });

  return (
    <Dialog fullScreen={fullScreen} open={isOpen} onClose={onClose}>
      <DialogTitle>
        {jobId && hasJob ? `Logs of ${jobs.items[jobId].name}` : 'Logs'}
      </DialogTitle>
      <DialogContent>
        {jobId && hasJob ? (
          <>
            <pre>{jobs.items[jobId].log}</pre>
            <div ref={logRef} />
          </>
        ) : (
          <div style={{ textAlign: 'center' }}>
            <CircularProgress />
          </div>
        )}
      </DialogContent>
      <DialogActions>
        {timer && (
          <>
            <Icon className="fas fa-sync fa-spin" />
            <TextField
              label="Refresh every"
              variant="filled"
              value={timeout}
              onChange={handleTimeoutChange}
              onBlur={handleTimeoutChange}
              InputProps={{
                endAdornment: <InputAdornment position="end">s</InputAdornment>
              }}
            />
          </>
        )}
        <Button onClick={internalOnClose} color="primary" autoFocus>
          Close
        </Button>
      </DialogActions>
    </Dialog>
  );
}

// $FlowFixMe: flow disabled for this line
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
