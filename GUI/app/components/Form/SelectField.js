// @flow
import React from 'react';
import { useField } from 'formik';
import { makeStyles } from '@material-ui/core/styles';
import InputLabel from '@material-ui/core/InputLabel';
import MenuItem from '@material-ui/core/MenuItem';
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import type { SimpleMapType } from '../../types/common';

const useStyles = makeStyles(theme => ({
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120
  }
}));

const OPTION_MAPPER = ([k, v]) => (
  <MenuItem key={k} value={k}>
    {v}
  </MenuItem>
);
const EMPTY_OPTION = emptyText => (
  <MenuItem key="__EMPTY__" value="">
    {emptyText}
  </MenuItem>
);

type SelectFieldProp = {
  label: string,
  name: string,
  options: SimpleMapType<string>,
  helperText?: ?string,
  emptyText?: string,
  addEmpty?: boolean,
  required?: boolean,
  multiple?: boolean
};

SelectField.defaultProps = {
  emptyText: 'None',
  helperText: null,
  addEmpty: false,
  required: false,
  multiple: false
};

export default function SelectField({
  options,
  addEmpty,
  emptyText,
  helperText,
  label,
  required,
  ...props
}: SelectFieldProp) {
  const classes = useStyles();
  const [
    { name, onBlur, onChange, value, multiple },
    { error, touched }
    // $FlowFixMe
  ] = useField({ as: 'select', ...props });
  const entries = Object.entries(options).map(OPTION_MAPPER);
  if (addEmpty) {
    entries.unshift(EMPTY_OPTION(emptyText));
  }
  const hasHelperText = !!((touched && error) || helperText);
  const displayedText = touched && error ? error : helperText;
  return (
    <FormControl
      className={classes.formControl}
      fullWidth
      error={!!(touched && error)}
      required={required}
    >
      <InputLabel shrink>{label}</InputLabel>
      <Select
        name={name}
        value={value}
        onBlur={onBlur}
        onChange={onChange}
        displayEmpty={addEmpty}
        multiple={multiple}
      >
        {entries}
      </Select>
      {hasHelperText ? <FormHelperText>{displayedText}</FormHelperText> : null}
    </FormControl>
  );
}
