// @flow
import * as React from 'react';
import { useField } from 'formik';
import { makeStyles } from '@material-ui/core/styles';
import FormGroup from '@material-ui/core/FormGroup';
import SwitchField from './SwitchField';
import TextField from './TextField';
import type { MapType } from '../../types/common';

const useStyles = makeStyles(theme => ({
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120,
    flexGrow: 1
  },
  formIndented: {
    margin: theme.spacing(1),
    marginTop: theme.spacing(-2),
    flexGrow: 1
  }
}));

type CustomArgumentsFieldProps = {
  name: string,
  labelEnable: string,
  labelTextField?: string
};

CustomArgumentsField.defaultProps = {
  labelTextField: 'Custom Arguments'
};

export type CustomArgumentsValue = {
  enable: boolean,
  value: string
};

export function ProcessCustomArguments(
  params: *,
  customArguments: CustomArgumentsValue,
  parameterName: string
) {
  const tmp = { ...params };
  if (customArguments && customArguments.enable && customArguments.value) {
    tmp[parameterName] = customArguments.value;
  }
  return tmp;
}

export default function CustomArgumentsField({
  name,
  labelEnable,
  labelTextField
}: CustomArgumentsFieldProps) {
  const classes = useStyles();
  const [, meta] = useField(`${name}.enable`);
  const { value } = meta;
  return (
    <FormGroup row className={classes.formControl}>
      <>
        <SwitchField label={labelEnable} name={`${name}.enable`} />
        {!!value && (
          <FormGroup row className={classes.formIndented}>
            <TextField
              label={labelTextField || 'Custom Arguments'}
              name={`${name}.value`}
              type="text"
            />
          </FormGroup>
        )}
      </>
    </FormGroup>
  );
}
