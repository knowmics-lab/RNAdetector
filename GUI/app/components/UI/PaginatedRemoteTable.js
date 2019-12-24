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
import Snackbar from './Snackbar';
import type { Action, StateType } from '../../reducers/types';

export type TableColumn = {
  id: string,
  label: string,
  minWidth?: number,
  align?: 'left' | 'right' | 'center' | 'justify',
  format?: ({ +[string]: * }) => mixed
};

export type DispatchActions = {
  requestPage: number => Action,
  changeRowsPerPage: number => Action,
  resetErrorMessage: () => Action
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
  size?: 'small' | 'medium',
  columns: TableColumn[],
  keyField?: string,
  requestPage: number => void,
  changeRowsPerPage: number => void,
  resetErrorMessage: () => void,
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

type TableState = {
  showError: boolean
};

class PaginatedRemoteTable extends Component<TableProps, TableState> {
  static defaultProps = {
    keyField: 'id',
    size: 'small'
  };

  constructor(props) {
    super(props);
    this.state = {
      showError: false
    };
  }

  static getDerivedStateFromProps(props: TableProps, state: TableState) {
    const { showError: oldShowError } = state;
    const { paginationState } = props;
    if (!oldShowError && paginationState.error) {
      return {
        ...state,
        showError: true
      };
    }
    return null;
  }

  handleChangePage = (event, newPage) => {
    const { requestPage } = this.props;
    requestPage(newPage + 1);
  };

  handleChangeRowsPerPage = event => {
    const { changeRowsPerPage } = this.props;
    changeRowsPerPage(+event.target.value);
  };

  handleClose = () => {
    const { resetErrorMessage } = this.props;
    resetErrorMessage();
    this.setState({
      showError: false
    });
  };

  render() {
    const { showError } = this.state;

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
      fetching,
      fetched
    } = paginationState;

    if (
      !currentPage &&
      !totalRows &&
      !fetching &&
      !fetched &&
      !paginationState.error
    ) {
      requestPage(1);
    }

    const isLoading = fetching || !totalRows;
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
        <Snackbar
          message={`An error occurred: ${paginationState.error || ''}!`}
          isOpen={showError}
          setClosed={this.handleClose}
          variant="error"
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
