// @flow
import React from 'react';
import has from 'lodash/has';
import Paper from '@material-ui/core/Paper';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableHead from '@material-ui/core/TableHead';
import TablePagination from '@material-ui/core/TablePagination';
import TableRow from '@material-ui/core/TableRow';
import TableContainer from '@material-ui/core/TableContainer';
import LinearProgress from '@material-ui/core/LinearProgress';
import { makeStyles } from '@material-ui/core/styles';
import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import type { ComponentType } from 'react';
import type { StatePaginationType } from '../../types/common';
import Snackbar from './Snackbar';
import type { Action, StateType } from '../../reducers/types';

export type TableColumn = {
  id: string,
  label: string,
  minWidth?: number,
  align?: 'left' | 'right' | 'center' | 'justify',
  format?: (*) => mixed
};

export type DispatchActions = {
  requestPage: number => Action,
  changeRowsPerPage: number => Action,
  handleSnackbarClose: () => Action
};

export type StateProps = {
  paginationState: StatePaginationType,
  pagesCollection: { +[number]: { +[string]: * }[] }
};

export type ConnectedTableProps = {
  columns: TableColumn[],
  keyField?: string
};

export type TableProps = {
  columns: TableColumn[],
  keyField?: string,
  requestPage: number => void,
  changeRowsPerPage: number => void,
  handleSnackbarClose: () => void
} & StateProps;

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%'
  },
  container: {
    maxHeight: 440
  },
  stickyStyle: {
    backgroundColor: 'white'
  },
  loading: {
    width: '100%',
    '& > * + *': {
      marginTop: theme.spacing(2)
    }
  }
}));

PaginatedRemoteTable.defaultProps = {
  keyField: 'id'
};

export default function PaginatedRemoteTable({
  columns,
  keyField,
  paginationState,
  pagesCollection,
  requestPage,
  changeRowsPerPage,
  handleSnackbarClose
}: TableProps) {
  const classes = useStyles();

  const handleChangePage = (event, newPage) => {
    requestPage(newPage + 1);
  };
  const handleChangeRowsPerPage = event => {
    changeRowsPerPage(+event.target.value);
  };
  const handleClose = () => {
    handleSnackbarClose();
  };

  const {
    current_page: currentPage,
    per_page: rowsPerPage,
    total: totalRows,
    isFetching,
    isError,
    errorMessage
  } = paginationState;

  if (!currentPage && !totalRows && !isFetching && !isError) {
    requestPage(1);
  }

  const isLoading = isFetching || !totalRows;
  const data =
    !isLoading && has(pagesCollection, currentPage)
      ? pagesCollection[currentPage]
      : [];

  return (
    <Paper className={classes.root}>
      {isLoading && (
        <div className={classes.loading}>
          <LinearProgress />
        </div>
      )}
      <TableContainer className={classes.container}>
        <Table stickyHeader size="medium">
          <TableHead>
            <TableRow>
              {columns.map(column => (
                <TableCell
                  key={column.id}
                  align={column.align}
                  style={{ minWidth: column.minWidth }}
                  className={classes.stickyStyle}
                >
                  {column.label}
                </TableCell>
              ))}
            </TableRow>
          </TableHead>
          <TableBody>
            {!isLoading &&
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
              ))}
          </TableBody>
        </Table>
      </TableContainer>
      <TablePagination
        rowsPerPageOptions={[1, 15, 30, 50, 100]}
        component="div"
        count={totalRows || 0}
        rowsPerPage={rowsPerPage || 15}
        page={(currentPage || 1) - 1}
        onChangePage={handleChangePage}
        onChangeRowsPerPage={handleChangeRowsPerPage}
      />
      <Snackbar
        message={`An error occurred: ${errorMessage}!`}
        isOpen={isError}
        setClosed={handleClose}
        variant="error"
      />
    </Paper>
  );
}

export function ConnectTable(
  mapStateToProps: StateType => StateProps,
  dispatchProps: DispatchActions
): ComponentType<ConnectedTableProps> {
  const mapDispatchToProps = dispatch =>
    bindActionCreators(dispatchProps, dispatch);

  // $FlowFixMe: flow disabled for this line
  return connect(mapStateToProps, mapDispatchToProps)(PaginatedRemoteTable);
}
