// @flow
import React from 'react';
import clsx from 'clsx';
import CheckCircleIcon from '@material-ui/icons/CheckCircle';
import ErrorIcon from '@material-ui/icons/Error';
import InfoIcon from '@material-ui/icons/Info';
import CloseIcon from '@material-ui/icons/Close';
import { amber, green } from '@material-ui/core/colors';
import IconButton from '@material-ui/core/IconButton';
import { Snackbar as MaterialSnackbar } from '@material-ui/core';
import SnackbarContent from '@material-ui/core/SnackbarContent';
import WarningIcon from '@material-ui/icons/Warning';
import { makeStyles } from '@material-ui/core/styles';

const variantIcon = {
  success: CheckCircleIcon,
  warning: WarningIcon,
  error: ErrorIcon,
  info: InfoIcon
};

const useStylesWrapper = makeStyles(theme => ({
  success: {
    backgroundColor: green[600]
  },
  error: {
    backgroundColor: theme.palette.error.dark
  },
  info: {
    backgroundColor: theme.palette.primary.main
  },
  warning: {
    backgroundColor: amber[700]
  },
  icon: {
    fontSize: 20
  },
  iconVariant: {
    opacity: 0.9,
    marginRight: theme.spacing(1)
  },
  message: {
    display: 'flex',
    alignItems: 'center'
  }
}));

type ContentWrapperProps = {
  message: string,
  onClose: (*, string) => void,
  variant: 'success' | 'warning' | 'error' | 'info'
};

function ContentWrapper({ message, onClose, variant }: ContentWrapperProps) {
  const classes = useStylesWrapper();
  const Icon = variantIcon[variant];

  return (
    <SnackbarContent
      className={clsx(classes[variant])}
      message={
        <span className={classes.message}>
          <Icon className={clsx(classes.icon, classes.iconVariant)} />
          {message}
        </span>
      }
      action={[
        <IconButton
          key="close"
          aria-label="close"
          color="inherit"
          onClick={onClose}
        >
          <CloseIcon className={classes.icon} />
        </IconButton>
      ]}
    />
  );
}

export type SnackbarProps = {
  message: string,
  isOpen: boolean,
  setClosed: () => void,
  variant: 'success' | 'warning' | 'error' | 'info',
  duration?: number,
  anchorVertical?: 'top' | 'bottom',
  anchorHorizontal?: 'left' | 'center' | 'right'
};

Snackbar.defaultProps = {
  duration: 3000,
  anchorVertical: 'bottom',
  anchorHorizontal: 'right'
};

export default function Snackbar({
  message,
  isOpen,
  setClosed,
  variant,
  duration,
  anchorVertical,
  anchorHorizontal
}: SnackbarProps) {
  const onClose = (event, reason) => {
    if (reason === 'clickaway') {
      return;
    }
    setClosed();
  };

  return (
    <MaterialSnackbar
      anchorOrigin={{
        vertical: anchorVertical,
        horizontal: anchorHorizontal
      }}
      open={isOpen}
      autoHideDuration={duration}
      onClose={onClose}
    >
      <ContentWrapper onClose={onClose} variant={variant} message={message} />
    </MaterialSnackbar>
  );
}
