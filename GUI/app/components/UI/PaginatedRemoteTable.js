// @flow
import React, { Component } from 'react';
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
import { withStyles } from '@material-ui/core/styles';
import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import type { ComponentType } from 'react';
import type { StatePaginationType } from '../../types/common';
import type { StateType } from '../../reducers/types';

export type TableColumn = {
  id: string,
  label: string,
  minWidth?: number,
  align?: 'left' | 'right' | 'center' | 'justify',
  format?: ({ +[string]: * }) => mixed
};

export type DispatchActions = {
  requestPage: number => *,
  changeRowsPerPage: number => *
};

export type StateProps = {
  paginationState: StatePaginationType,
  pagesCollection: { +[number]: { +[string]: * }[] }
};

export type ConnectedTableProps = {
  size?: 'small' | 'medium',
  columns: TableColumn[],
  keyField?: string,
  onPageChange: number => void
};

export type TableProps = {
  size?: 'small' | 'medium',
  columns: TableColumn[],
  keyField?: string,
  onPageChange: number => void,
  requestPage: number => void,
  changeRowsPerPage: number => void,
  paginationState: StatePaginationType,
  pagesCollection: { +[number]: { +[string]: * }[] },
  classes: {
    root: *,
    container: *,
    stickyStyle: *,
    loading: *
  }
};

const styles = theme => ({
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
});

class PaginatedRemoteTable extends Component<TableProps> {
  props: TableProps;

  static defaultProps = {
    keyField: 'id',
    size: 'small'
  };

  componentDidUpdate(prevProps: TableProps) {
    const {
      paginationState: { current_page: prevPage }
    } = prevProps;
    const {
      paginationState: { current_page: currentPage },
      onPageChange
    } = this.props;
    if (currentPage && currentPage !== prevPage) onPageChange(currentPage);
  }

  handleChangePage = (event, newPage) => {
    const { requestPage, onPageChange } = this.props;
    requestPage(newPage + 1);
    onPageChange(newPage + 1);
  };

  handleChangeRowsPerPage = event => {
    const { changeRowsPerPage } = this.props;
    changeRowsPerPage(+event.target.value);
  };

  render() {
    const {
      size,
      classes,
      columns,
      keyField,
      paginationState,
      pagesCollection,
      requestPage
    } = this.props;

    const {
      current_page: currentPage,
      per_page: rowsPerPage,
      total: totalRows,
      fetching
    } = paginationState;

    if (currentPage === null && totalRows === null && !fetching) {
      requestPage(1);
    }
    const isLoading = fetching || totalRows === null;
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
          <Table stickyHeader size={size}>
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
                          {column.format ? column.format(row) : value}
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
          onChangePage={this.handleChangePage}
          onChangeRowsPerPage={this.handleChangeRowsPerPage}
        />
      </Paper>
    );
  }
}

export default withStyles(styles)(PaginatedRemoteTable);

export function ConnectTable(
  mapStateToProps: StateType => StateProps,
  dispatchProps: DispatchActions
): ComponentType<ConnectedTableProps> {
  const mapDispatchToProps = dispatch =>
    bindActionCreators(dispatchProps, dispatch);

  // $FlowFixMe: flow disabled for this line
  return connect(
    mapStateToProps,
    mapDispatchToProps
  )(withStyles(styles)(PaginatedRemoteTable));
}
