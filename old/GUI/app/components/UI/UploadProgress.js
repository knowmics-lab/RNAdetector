// @flow
import * as React from 'react';
import { makeStyles } from '@material-ui/core/styles';
import Box from '@material-ui/core/Box';
import FormGroup from '@material-ui/core/FormGroup';
import Grid from '@material-ui/core/Grid';
import { LinearProgress } from '@material-ui/core';
import byteSize from 'byte-size';

const useStyles = makeStyles(theme => ({
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120
  }
}));

type Props = {
  isUploading: boolean,
  uploadFile: string,
  uploadedBytes: number,
  uploadedPercent: number,
  uploadTotal: number
};

export default function UploadProgress({
  isUploading,
  uploadFile,
  uploadedBytes,
  uploadedPercent,
  uploadTotal
}: Props) {
  const classes = useStyles();
  return (
    <>
      {isUploading && (
        <FormGroup row className={classes.formControl}>
          <Grid container justify="center" alignItems="center" spacing={1}>
            <Grid item xs={12}>
              <Box fontWeight="fontWeightMedium">{`Uploading ${uploadFile}...`}</Box>
            </Grid>
            <Grid item xs={12}>
              <Grid
                container
                justify="space-evenly"
                alignItems="center"
                spacing={1}
              >
                <Grid item xs={10}>
                  <LinearProgress
                    variant="determinate"
                    value={uploadedPercent}
                  />
                </Grid>
                <Grid item xs>
                  {`${byteSize(uploadedBytes)}/${byteSize(uploadTotal)}`}
                </Grid>
              </Grid>
            </Grid>
          </Grid>
        </FormGroup>
      )}
    </>
  );
}
