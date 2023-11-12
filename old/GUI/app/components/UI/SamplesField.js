/* eslint-disable no-restricted-syntax */
// @flow
import * as React from 'react';
import path from 'path';
import { activeWindow, api } from 'electron-util';
import { makeStyles } from '@material-ui/core/styles';
import ListItemSecondaryAction from '@material-ui/core/ListItemSecondaryAction';
import IconButton from '@material-ui/core/IconButton';
import ListSubheader from '@material-ui/core/ListSubheader';
import Icon from '@material-ui/core/Icon';
import FileType from 'file-type';
import { FieldArray, useField, useFormikContext } from 'formik';
import Grid from '@material-ui/core/Grid';
import { Hidden, InputLabel } from '@material-ui/core';
import Typography from '@material-ui/core/Typography';
import Button from '@material-ui/core/Button';
import Box from '@material-ui/core/Box';
import FormGroup from '@material-ui/core/FormGroup';
import * as Api from '../../api';
import type { File } from './FileSelector2';
import type { FileFilter } from '../../types/common';
import FileSelector from './FileSelector2';
import TextField from '../Form/TextField';
import MultiDropZone from './MultiDropZone';

const EMPTY_SAMPLE = { code: '', first: undefined, second: undefined };

const useStyles = makeStyles(theme => ({
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120
  }
}));

type SampleData = {
  code: string,
  first: ?(File[]),
  second: ?(File[])
};

type Props = {
  name: string,
  code: string,
  inputType: string,
  paired: boolean,
  disabled: boolean
};

function isEmptyFile(file: ?(File[])): boolean %checks {
  return !file || (file && !file[0]);
}

function isEmptySample(sample: SampleData): boolean %checks {
  return (
    !sample.code && isEmptyFile(sample.first) && isEmptyFile(sample.second)
  );
}

function fileExists(file: File, samples: SampleData[]) {
  for (const sample of samples) {
    if (!isEmptySample(sample)) {
      if (!isEmptyFile(sample.first) && sample.first[0].path === file.path)
        return true;
      if (!isEmptyFile(sample.second) && sample.second[0].path === file.path)
        return true;
    }
  }
  return false;
}

function getCode(file: File, paired: boolean): [?string, boolean] {
  const basename = file.name.split('.')[0];
  const p = basename.indexOf('_');
  if (paired && p < 0) return [undefined, true];
  const code = paired ? basename.substring(0, p) : basename;
  const first = paired
    ? basename
        .substring(p)
        .toUpperCase()
        .includes('_R1')
    : true;
  return [code, first];
}

export function prepareSamplesArray(samples: SampleData[]): { code: string }[] {
  return samples.map(s => {
    return {
      code: s.code
    };
  });
}

export function prepareFileArray(samples: SampleData[]): [?File, ?File][] {
  return samples.map(s => {
    if (isEmptySample(s)) return [null, null];
    return [
      isEmptyFile(s.first) ? null : s.first[0],
      isEmptyFile(s.second) ? null : s.second[0]
    ];
  });
}

SamplesField.defaultProps = {};

export default function SamplesField({
  name,
  code,
  inputType,
  paired,
  disabled
}: Props) {
  const classes = useStyles();
  const [, meta, fieldHelpers] = useField<SampleData[]>(name);
  const { value: samples } = meta;
  const { setValue } = fieldHelpers;
  const single = samples.length <= 1;

  const removeHandler = React.useCallback(
    (helpers, i) => e => {
      e.preventDefault();
      helpers.remove(i);
    },
    []
  );

  const makeSampleForm = React.useCallback(
    (helpers, i) => {
      return (
        <FormGroup row className={classes.formControl} key={`sample-${i}`}>
          <Grid
            container
            justify="space-around"
            alignItems="center"
            spacing={3}
          >
            {!single && (
              <Grid item md={3}>
                <TextField
                  label="Sample Code"
                  name={`${name}.${i}.code`}
                  placeholder={`${code}_${i + 1}`}
                />
              </Grid>
            )}
            {single && (
              <Grid item md={3}>
                <Typography
                  variant="h6"
                  align="right"
                  style={{ fontSize: '1rem' }}
                >
                  {`Choose ${inputType.toUpperCase()} file:`}
                </Typography>
              </Grid>
            )}
            <Grid item md={paired ? 4 : 8}>
              <FileSelector
                name={`${name}.${i}.first`}
                filters={Api.Utils.analysisFileExtensions(inputType)}
                disabled={disabled}
              />
            </Grid>
            {paired && (
              <Grid item md={4}>
                <FileSelector
                  name={`${name}.${i}.second`}
                  filters={Api.Utils.analysisFileExtensions(inputType)}
                  disabled={disabled}
                />
              </Grid>
            )}
            {!single && (
              <Grid item xs={1}>
                <IconButton onClick={removeHandler(helpers, i)}>
                  <Icon className="fas fa-trash" />
                </IconButton>
              </Grid>
            )}
          </Grid>
        </FormGroup>
      );
    },
    [name, single, classes, paired, disabled, removeHandler]
  );

  const addHandler = React.useCallback(
    helpers => (e: MouseEvent) => {
      e.preventDefault();
      helpers.push({ ...EMPTY_SAMPLE });
    },
    []
  );

  const bulkAddFileHandler = React.useCallback(
    () => (newFiles: File[]) => {
      const cleanSamples = samples.filter(s => !isEmptySample(s));
      const insertUpdateOrOverwrite = (
        sampleCode: string,
        newSample: $Shape<SampleData>
      ) => {
        for (let i = 0; i < cleanSamples.length; i += 1) {
          if (cleanSamples[i].code === sampleCode) {
            cleanSamples[i] = {
              ...cleanSamples[i],
              ...newSample
            };
            return;
          }
        }
        cleanSamples.push({
          ...EMPTY_SAMPLE,
          ...newSample
        });
      };

      const toAddFiles = newFiles.filter(f => !fileExists(f, cleanSamples));
      if (toAddFiles.length > 0) {
        toAddFiles.forEach(f => {
          const [sCode, first] = getCode(f, paired);
          if (sCode) {
            if (first) {
              insertUpdateOrOverwrite(sCode, {
                code: sCode,
                first: [f]
              });
            } else {
              insertUpdateOrOverwrite(sCode, {
                code: sCode,
                second: [f]
              });
            }
          }
        });
        setValue(cleanSamples);
      }
    },
    [samples, paired]
  );

  return (
    <FieldArray
      name="samples"
      render={helpers => (
        <>
          {samples.map((s, i) => makeSampleForm(helpers, i))}
          <FormGroup row className={classes.formControl}>
            <MultiDropZone
              onFileAdd={bulkAddFileHandler()}
              onSampleAdd={addHandler(helpers)}
              filters={Api.Utils.analysisFileExtensions(inputType)}
            />
          </FormGroup>
        </>
      )}
    />
  );
}
