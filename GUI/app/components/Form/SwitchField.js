// @flow
import React from 'react';
import { useField } from 'formik';
import { makeStyles } from '@material-ui/core/styles';
import FormControl from '@material-ui/core/FormControl';
import FormGroup from '@material-ui/core/FormGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormHelperText from '@material-ui/core/FormHelperText';
import Switch from '@material-ui/core/Switch';

const useStyles = makeStyles(theme => ({
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120
  }
}));

export type SwitchFieldProps = {
  label: string,
  name: string,
  required?: boolean
};

SwitchField.defaultProps = {
  required: false
};

export default function SwitchField({
  label,
  required,
  ...props
}: SwitchFieldProps) {
  const classes = useStyles();
  const [
    { name, checked, onBlur, onChange, value },
    { error, touched }
  ] = useField({ ...props, type: 'checkbox' });

  return (
    <FormControl
      className={classes.formControl}
      fullWidth
      error={!!(touched && error)}
    >
      <FormGroup>
        <FormControlLabel
          control={
            <Switch
              name={name}
              checked={checked}
              onChange={onChange}
              onBlur={onBlur}
              value={value}
            />
          }
          label={label}
        />
      </FormGroup>
      {touched && error ? <FormHelperText>{error}</FormHelperText> : null}
    </FormControl>
  );
}
