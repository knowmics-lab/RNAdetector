/* eslint-disable no-restricted-syntax */
// @flow
import * as React from 'react';
import path from 'path';
import { activeWindow, api } from 'electron-util';
import { makeStyles } from '@material-ui/core/styles';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemAvatar from '@material-ui/core/ListItemAvatar';
import ListItemSecondaryAction from '@material-ui/core/ListItemSecondaryAction';
import ListItemText from '@material-ui/core/ListItemText';
import Avatar from '@material-ui/core/Avatar';
import IconButton from '@material-ui/core/IconButton';
import ListSubheader from '@material-ui/core/ListSubheader';
import Icon from '@material-ui/core/Icon';
import mimetype2fa from 'mimetype-to-fontawesome';
import FileType from 'file-type';
import Box from '@material-ui/core/Box';
import Divider from '@material-ui/core/Divider';
import type { FileFilter } from '../../types/common';
import * as Api from '../../api';

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    maxWidth: 250,
    backgroundColor: theme.palette.background.paper
  },
  title: {
    whiteSpace: 'nowrap',
    overflow: 'hidden',
    textOverflow: 'ellipsis'
  }
}));

export type File = {
  name: string,
  path: string,
  type: string,
  ext: string
};

type Props = {
  title?: ?string,
  onFileAdd: (File[]) => void,
  onFileRemove: File => void,
  multiple?: boolean,
  disabled?: boolean,
  filters?: FileFilter[]
};

FileSelector.defaultProps = {
  title: null,
  multiple: false,
  disabled: false,
  filters: [{ name: 'All Files', extensions: ['*'] }]
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

export default function FileSelector({
  title,
  multiple,
  disabled,
  filters,
  onFileAdd,
  onFileRemove
}: Props) {
  const classes = useStyles();
  const [files, setFiles] = React.useState(([]: File[]));
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

  const fileExists = (file: File) => {
    for (let i = 0, l = files.length; i < l; i += 1) {
      if (files[i].path === file.path) return true;
    }
    return false;
  };
  const realAddFile = async filePaths => {
    if (filePaths.length === 0) return;
    const selectedFiles = (await processFiles(filePaths)).filter(
      f => isExtensionAllowed(f.ext, processedFilters) && !fileExists(f)
    );
    if (selectedFiles.length === 0) return;
    setFiles(prevFiles => [...prevFiles, ...selectedFiles]);
    onFileAdd(selectedFiles);
  };

  const handleAdd = async () => {
    const properties = [];
    if (multiple) properties.push('multiSelections');
    const { canceled, filePaths } = await api.dialog.showOpenDialog(
      activeWindow(),
      {
        filters,
        properties
      }
    );
    if (!canceled) {
      if (filePaths) {
        await realAddFile(filePaths);
      }
    }
  };

  const handleRemove = f => () => {
    setFiles(prevFiles => prevFiles.filter(o => o.path !== f.path));
    onFileRemove(f);
  };

  const handlePrevent = e => {
    e.preventDefault();
  };

  const handleDrop = e => {
    e.preventDefault();
    if (!multiple && files.length < 1) {
      const fileArray = Api.Utils.toArray(e.dataTransfer.files);
      let filePaths = fileArray.map(t => t.path);
      if (!multiple) {
        filePaths = [filePaths.shift()];
      }
      realAddFile(filePaths)
        .then(() => true)
        .catch(() => false);
    }
    return false;
  };

  const getTitle = () => {
    if (title !== null) return title;
    return multiple ? 'Select Files' : 'Select a file';
  };

  const makeSubheader = () => {
    if (!multiple && files.length > 0) return null;
    return (
      <ListSubheader>
        {getTitle()}
        <ListItemSecondaryAction>
          {(multiple || (!multiple && files.length < 1)) && (
            <IconButton edge="end" onClick={handleAdd} disabled={disabled}>
              <Icon className="fas fa-plus" />
            </IconButton>
          )}
        </ListItemSecondaryAction>
      </ListSubheader>
    );
  };

  const m2f = mimetype2fa({ prefix: 'fa-' });
  return (
    <Box
      className={classes.root}
      boxShadow={1}
      borderRadius="borderRadius"
      onDragOver={handlePrevent}
      onDragEnter={handlePrevent}
      onDragLeave={handlePrevent}
      onDragEnd={handlePrevent}
      onDrop={handleDrop}
    >
      <List subheader={makeSubheader()} dense>
        {multiple && files.length > 0 && <Divider />}
        {files.map(f => (
          <ListItem key={f.path}>
            <ListItemAvatar>
              <Avatar>
                <i className={`fas ${m2f(f.type)}`} />
              </Avatar>
            </ListItemAvatar>
            <ListItemText
              className={classes.title}
              primary={f.name}
              title={f.name}
            />
            {!disabled && (
              <ListItemSecondaryAction>
                <IconButton
                  edge="end"
                  color="secondary"
                  onClick={handleRemove(f)}
                >
                  <i className="fas fa-times" />
                </IconButton>
              </ListItemSecondaryAction>
            )}
          </ListItem>
        ))}
      </List>
    </Box>
  );
}
