// @flow
import React from 'react';
import ReactDOM from 'react-dom';
import { useField, useFormikContext } from 'formik';
import { makeStyles } from '@material-ui/core/styles';
import { TextField as MaterialTextField, Icon } from '@material-ui/core';
import { remote } from 'electron';
import IconButton from '@material-ui/core/IconButton';
import Input from '@material-ui/core/Input';
import InputLabel from '@material-ui/core/InputLabel';
import InputAdornment from '@material-ui/core/InputAdornment';
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';

const useStyles = makeStyles(theme => ({
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120
  }
}));

type TextFieldProps = {
  label: string,
  name: string,
  required?: boolean
};

FileField.defaultProps = {
  required: false
};

export default function FileField({
  label,
  required,
  ...props
}: TextFieldProps) {
  const classes = useStyles();
  const { setFieldValue } = useFormikContext();
  const [{ name, onBlur, onChange, value }, { error, touched }] = useField(
    props
  );
  return (
    <FormControl
      className={classes.formControl}
      fullWidth
      error={!!(touched && error)}
      required={required}
    >
      <InputLabel>{label}</InputLabel>
      <Input
        name={name}
        type="text"
        value={value}
        onChange={onChange}
        onBlur={onBlur}
        endAdornment={
          <InputAdornment position="end">
            <IconButton
            // onClick={handleClickShowPassword}
            // onMouseDown={handleMouseDownPassword}
            >
              <Icon className="fas fa-ellipsis-h" />
            </IconButton>
          </InputAdornment>
        }
      />
      {touched && error ? <FormHelperText>{error}</FormHelperText> : null}
    </FormControl>
  );
}
