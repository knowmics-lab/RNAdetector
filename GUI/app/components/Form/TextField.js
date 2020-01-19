// @flow
import * as React from 'react';
import { useField } from 'formik';
import { makeStyles } from '@material-ui/core/styles';
import { TextField as MaterialTextField } from '@material-ui/core';

const useStyles = makeStyles(theme => ({
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120
  }
}));

type TextFieldProps = {
  label: string,
  name: string,
  type?: string,
  helperText?: React.Node,
  placeholder?: string,
  readOnly?: boolean,
  required?: boolean,
  multiline?: boolean
};

TextField.defaultProps = {
  type: 'text',
  placeholder: '',
  helperText: null,
  required: false,
  readOnly: false,
  multiline: false
};

export default function TextField({
  label,
  placeholder,
  required,
  helperText,
  ...props
}: TextFieldProps) {
  const classes = useStyles();
  const [{ name, onBlur, onChange, value }, { error, touched }] = useField(
    props
  );
  return (
    <MaterialTextField
      fullWidth
      className={classes.formControl}
      label={label}
      required={required}
      error={!!(touched && error)}
      placeholder={placeholder}
      name={name}
      onBlur={onBlur}
      onChange={onChange}
      value={value}
      helperText={touched && error ? error : helperText}
    />
  );
}
