/* eslint-disable no-restricted-syntax */
// @flow
import * as React from 'react';
import path from 'path';
import { makeStyles } from '@material-ui/core/styles';
import Icon from '@material-ui/core/Icon';
import FileType from 'file-type';
import Grid from '@material-ui/core/Grid';
import { Hidden } from '@material-ui/core';
import Typography from '@material-ui/core/Typography';
import Button from '@material-ui/core/Button';
import * as Api from '../../api';
import type { File } from './FileSelector2';
import type { FileFilter } from '../../types/common';

const useStyles = makeStyles(theme => ({
  dropZoneContainer: {
    margin: theme.spacing(2)
  },
  gridContent: {
    padding: theme.spacing(2),
    textAlign: 'center'
  },
  dropZone: {
    textAlign: 'center',
    width: '100%',
    border: 'dashed',
    cursor: 'grab',
    overflow: 'hidden',
    position: 'relative',
    boxSizing: 'border-box',
    minHeight: '50px',
    borderColor: 'rgba(0, 0, 0, 0.12)',
    borderRadius: '4px'
  },
  'dropZone:active': {
    cursor: 'grabbing'
  }
}));

type Props = {
  onSampleAdd: MouseEvent => void,
  onFileAdd: (File[]) => void,
  filters?: FileFilter[]
};

function processFiltersList(filters: ?(FileFilter[])): string[] {
  const res = [];
  for (const filter of filters || []) {
    for (const ext of filter.extensions) {
      if (ext === '*') res.unshift('*');
      else res.push(`.${ext.toLowerCase()}`);
    }
  }
  return res;
}

function isExtensionAllowed(ext: string, filters: string[]): boolean {
  if (filters.length && filters[0] === '*') return true;
  const lowerExt = ext.toLowerCase();
  return filters.includes(lowerExt);
}

MultiDropZone.defaultProps = {
  filters: [{ name: 'All Files', extensions: ['*'] }]
};

export default function MultiDropZone({
  onSampleAdd,
  onFileAdd,
  filters
}: Props) {
  const classes = useStyles();
  const processedFilters = React.useMemo(() => processFiltersList(filters), [
    filters
  ]);

  const processFiles = async (paths: string[]): Promise<File[]> => {
    return Promise.all(
      paths.map(async f => {
        const t = await FileType.fromFile(f);
        return {
          name: path.basename(f),
          path: f,
          type: (t && t.mime) || 'text/plain',
          ext: path.extname(f)
        };
      })
    );
  };

  const doAddFiles = async filePaths => {
    if (filePaths.length === 0) return;
    const selectedFiles = (await processFiles(filePaths)).filter(f =>
      isExtensionAllowed(f.ext, processedFilters)
    );
    if (selectedFiles.length === 0) return;
    onFileAdd(selectedFiles);
  };

  const handlePrevent = e => {
    e.preventDefault();
  };

  const handleDrop = e => {
    e.preventDefault();
    const fileArray = Api.Utils.toArray(e.dataTransfer.files);
    doAddFiles(fileArray.map(t => t.path))
      .then(() => true)
      .catch(() => false);
    return false;
  };

  return (
    <Grid
      container
      direction="row"
      alignItems="center"
      spacing={3}
      className={classes.dropZoneContainer}
    >
      <Hidden only="sm">
        <Grid
          item
          md={5}
          zeroMinWidth
          className={classes.dropZone}
          onDragOver={handlePrevent}
          onDragEnter={handlePrevent}
          onDragLeave={handlePrevent}
          onDragEnd={handlePrevent}
          onDrop={handleDrop}
        >
          <Typography variant="h6">
            <Icon className="fas fa-upload" />
            &nbsp; Drag and drop files here
          </Typography>
        </Grid>
        <Grid item md={2} className={classes.gridContent}>
          <Typography variant="h6">or</Typography>
        </Grid>
      </Hidden>
      <Grid item xs={12} md={5} className={classes.gridContent}>
        <Button variant="outlined" onClick={onSampleAdd}>
          <Icon className="fas fa-plus" /> Add a sample manually
        </Button>
      </Grid>
    </Grid>
  );
}
