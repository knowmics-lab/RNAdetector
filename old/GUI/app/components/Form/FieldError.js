// @flow
import * as React from 'react';
import { useField } from 'formik';
import { makeStyles } from '@material-ui/core/styles';
import FormControl from '@material-ui/core/FormControl';
import FormLabel from '@material-ui/core/FormLabel';
import FormHelperText from '@material-ui/core/FormHelperText';

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    margin: theme.spacing(1),
    minWidth: 120
  },
  label: {
    fontSize: '0.75rem'
  }
}));

type CustomFieldProps = {
  name: string,
  label: string,
  children: React.ChildrenArray<*> | (number => React.Node),
  required?: boolean,
  helperText?: React.Node
};

CustomField.defaultProps = {
  required: false,
  helperText: null
};

export default function CustomField({
  label,
  required,
  children,
  helperText,
  ...props
}: CustomFieldProps) {
  const classes = useStyles();
  const [, { error, touched }] = useField(props);
  const hasError = !!(touched && error && typeof error === 'string');
  const finalHelperText = hasError ? error : helperText;
  return (
    <FormControl
      fullWidth
      required={required}
      className={classes.root}
      error={hasError}
    >
      <FormLabel className={classes.label}>{label}</FormLabel>
      {children}
      {!!finalHelperText && <FormHelperText>{finalHelperText}</FormHelperText>}
    </FormControl>
  );
}
