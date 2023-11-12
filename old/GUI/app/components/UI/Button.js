// @flow
import React, { useMemo, forwardRef } from 'react';
import { Link as RouterLink } from 'react-router-dom';
import { Button as OB } from '@material-ui/core';
import CircularProgress from '@material-ui/core/CircularProgress';
import makeStyles from '@material-ui/core/styles/makeStyles';
import { green } from '@material-ui/core/colors';

export type ButtonType = {
  color?: 'default' | 'inherit' | 'primary' | 'secondary',
  children: *,
  disabled?: boolean,
  href?: ?string,
  size?: 'small' | 'medium' | 'large',
  variant?: 'text' | 'outlined' | 'contained',
  onClick?: ?() => void
};

export default function Button({
  color: co,
  children: c,
  disabled: d,
  href,
  size: s,
  variant: v,
  onClick
}: ButtonType) {
  if (!href && onClick) {
    return (
      <OB variant={v} color={co} disabled={d} onClick={onClick} size={s}>
        {c}
      </OB>
    );
  }
  if (!onClick && href) {
    const renderLink = useMemo(
      () =>
        forwardRef((itemProps, ref) => (
          <RouterLink
            to={href}
            // eslint-disable-next-line react/jsx-props-no-spreading
            {...itemProps}
            innerRef={ref}
          />
        )),
      [href]
    );
    return (
      <OB variant={v} color={co} disabled={d} size={s} component={renderLink}>
        {c}
      </OB>
    );
  }
  throw new Error('You must specify onClick or href.');
}

Button.defaultProps = {
  color: 'default',
  size: 'medium',
  variant: 'text',
  disabled: false,
  href: null,
  onClick: null
};

const useStyles = makeStyles(theme => ({
  buttonWrapper: {
    margin: theme.spacing(1),
    position: 'relative'
  },
  buttonProgress: {
    color: green[500],
    position: 'absolute',
    top: '50%',
    left: '50%',
    marginTop: -12,
    marginLeft: -12
  }
}));

export type SubmitButtonType = {
  isSaving: boolean,
  text?: ?string
};

export function SubmitButton({ isSaving, text }: SubmitButtonType) {
  const classes = useStyles();
  return (
    <div className={classes.buttonWrapper}>
      <OB type="submit" variant="contained" color="primary" disabled={isSaving}>
        {text || 'Save'}
      </OB>
      {isSaving && (
        <CircularProgress size={24} className={classes.buttonProgress} />
      )}
    </div>
  );
}

SubmitButton.defaultProps = {
  text: 'Save'
};
