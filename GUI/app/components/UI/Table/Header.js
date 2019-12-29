// @flow
import React from 'react';
import { makeStyles } from '@material-ui/core/styles';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import TableCell from '@material-ui/core/TableCell';
import type { TableColumn } from './types';

const useStyles = makeStyles({
  stickyStyle: {
    backgroundColor: 'white'
  }
});

type Props = {
  columns: TableColumn[]
};

export default function Header({ columns }: Props) {
  const classes = useStyles();
  return (
    <TableHead>
      <TableRow>
        {columns.map(column =>
          column !== 'actions' ? (
            <TableCell
              key={column.id}
              align={column.align}
              style={{ minWidth: column.minWidth }}
              className={classes.stickyStyle}
            >
              {column.label}
            </TableCell>
          ) : (
            <TableCell
              key="actions"
              align="center"
              className={classes.stickyStyle}
            >
              Actions
            </TableCell>
          )
        )}
      </TableRow>
    </TableHead>
  );
}
