import React from 'react';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';
import useMediaQuery from '@material-ui/core/useMediaQuery';
import { useTheme } from '@material-ui/core/styles';
import { TextField } from '@material-ui/core';
import InputAdornment from '@material-ui/core/InputAdornment';

type LogsDialogProps = {
  jobId: number,
  open: boolean,
  onClose: () => void
};

function LogsDialog({ jobId, open, onClose }: LogsDialogProps) {
  const theme = useTheme();
  const [timeout, setTimeout] = React.useState(2);
  const fullScreen = useMediaQuery(theme.breakpoints.down('md'));

  const handleTimeoutChange = e => {
    console.log(e.target.value);
  };

  return (
    <>
      <Dialog fullScreen={fullScreen} open={open} onClose={onClose}>
        <DialogTitle>TODO</DialogTitle>
        <DialogContent>
          <DialogContentText>TODO</DialogContentText>
        </DialogContent>
        <DialogActions>
          <TextField
            label="Refresh every"
            variant="filled"
            value={timeout}
            onChange={handleTimeoutChange}
            endAdornment={<InputAdornment position="end">Kg</InputAdornment>}
          />
          <Button onClick={onClose} color="primary" autoFocus>
            Close
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
}

export default LogsDialog;
