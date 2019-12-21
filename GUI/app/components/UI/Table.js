// @flow
import React from 'react';
import {
  Table as MaterialTable,
  TableBody,
  TableCell,
  TableHead,
  TablePagination,
  TableRow,
  TableContainer,
  CircularProgress
} from '@material-ui/core';
import { makeStyles } from '@material-ui/core/styles';

export type TableColumn = {
  id: string,
  label: string,
  minWidth: number,
  align?: 'left' | 'right' | 'center' | 'justify',
  format?: mixed => string
};

export type TableProps = {
  columns: TableColumn[],
  keyField?: string,
  loading?: boolean,
  data: *,
  currentPage: number,
  totalRows: number,
  rowsPerPage: number,
  onPageChange: (page: number) => void,
  onRowsPerPageChange: (perPage: number) => void
};

const useStyles = makeStyles({
  root: {
    width: '100%'
  },
  container: {
    maxHeight: 440
  }
});

Table.defaultProps = {
  loading: false,
  keyField: 'id'
};

export default function Table({
  columns,
  keyField,
  loading,
  data,
  currentPage,
  totalRows,
  rowsPerPage,
  onPageChange,
  onRowsPerPageChange
}: TableProps) {
  const classes = useStyles();

  const handleChangePage = (event, newPage) => {
    onPageChange(newPage + 1);
  };

  const isLoading = loading || !data;

  const handleChangeRowsPerPage = event => {
    onRowsPerPageChange(+event.target.value);
  };

  return (
    <>
      <TableContainer className={classes.container}>
        <MaterialTable stickyHeader aria-label="sticky table">
          <TableHead>
            <TableRow>
              {columns.map(column => (
                <TableCell
                  key={column.id}
                  align={column.align}
                  style={{ minWidth: column.minWidth }}
                >
                  {column.label}
                </TableCell>
              ))}
            </TableRow>
          </TableHead>
          <TableBody>
            {isLoading ? (
              <TableRow>
                <TableCell align="center" colSpan={columns.length}>
                  <CircularProgress />
                </TableCell>
              </TableRow>
            ) : (
              data.map(row => (
                <TableRow
                  hover
                  role="checkbox"
                  tabIndex={-1}
                  key={row[keyField]}
                >
                  {columns.map(column => {
                    const value = row[column.id];
                    return (
                      <TableCell key={column.id} align={column.align}>
                        {column.format ? column.format(value) : value}
                      </TableCell>
                    );
                  })}
                </TableRow>
              ))
            )}
          </TableBody>
        </MaterialTable>
      </TableContainer>
      {!isLoading && (
        <TablePagination
          rowsPerPageOptions={[1, 10, 15, 20, 30, 100]}
          component="div"
          count={totalRows}
          rowsPerPage={rowsPerPage}
          page={currentPage - 1}
          onChangePage={handleChangePage}
          onChangeRowsPerPage={handleChangeRowsPerPage}
        />
      )}
    </>
  );
}
