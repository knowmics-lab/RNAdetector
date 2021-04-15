// @flow
import * as React from 'react';
import has from 'lodash/has';
import Paper from '@material-ui/core/Paper';
import Table from '@material-ui/core/Table';
import TablePagination from '@material-ui/core/TablePagination';
import TableContainer from '@material-ui/core/TableContainer';
import LinearProgress from '@material-ui/core/LinearProgress';
import { withStyles } from '@material-ui/core/styles';
import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import type { ComponentType } from 'react';
import type { SortingSpec, StatePaginationType } from '../../types/common';
import type { StateType } from '../../reducers/types';
import type {
  ReadOnlyData,
  RowActionType,
  TableColumn,
  ToolbarActionType
} from './Table/types';
import TableHeader from './Table/Header';
import TableBody from './Table/Body';
import TableToolbar from './Table/Toolbar';

export type DispatchActions = {
  requestPage: number => *,
  changeRowsPerPage: number => *,
  changeSorting: SortingSpec => *
};

export type StateProps = {
  paginationState: StatePaginationType,
  pagesCollection: { +[number]: { +[string]: * }[] }
};

export type ConnectedTableProps = {
  title?: ?(string | React.ChildrenArray<React.Element<*>>),
  size?: 'small' | 'medium',
  columns: TableColumn[],
  toolbar?: ToolbarActionType[],
  actions?: RowActionType[],
  keyField?: string,
  sortable?: boolean,
  onPageChange: number => void,
  autoRefresh?: boolean,
  autoRefreshCondition?: ?(ReadOnlyData[]) => boolean,
  autoRefreshAction?: ?(number) => void,
  autoRefreshTime?: number
};

export type TableProps = {
  title?: ?(string | React.ChildrenArray<React.Element<*>>),
  size?: 'small' | 'medium',
  columns: TableColumn[],
  toolbar?: ToolbarActionType[],
  actions?: RowActionType[],
  keyField?: string,
  sortable?: boolean,
  onPageChange: number => void,
  requestPage: number => void,
  changeRowsPerPage: number => void,
  autoRefresh?: boolean,
  autoRefreshCondition?: ?(ReadOnlyData[]) => boolean,
  autoRefreshAction?: ?(number) => void,
  autoRefreshTime?: number,
  changeSorting: SortingSpec => void,
  paginationState: StatePaginationType,
  pagesCollection: { +[number]: { +[string]: * }[] },
  classes: {
    root: *,
    container: *,
    stickyStyle: *,
    loading: *
  }
};

type TableState = {
  timer: *
};

const styles = theme => ({
  root: {
    width: '100%'
  },
  container: {
    // maxHeight: 440
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

class PaginatedRemoteTable extends React.Component<TableProps, TableState> {
  static defaultProps = {
    title: null,
    keyField: 'id',
    size: 'small',
    sortable: true,
    toolbar: [],
    actions: [],
    autoRefresh: false,
    autoRefreshCondition: undefined,
    autoRefreshAction: undefined,
    autoRefreshTime: 30000
  };

  constructor(props) {
    super(props);
    this.state = {
      timer: undefined
    };
  }

  componentDidMount(): void {
    const {
      paginationState: { current_page: currentPage },
      onPageChange,
      autoRefresh
    } = this.props;
    const { timer } = this.state;
    if (currentPage) onPageChange(currentPage);
    if (!timer && autoRefresh) {
      this.checkForRefresh();
    }
  }

  // noinspection JSCheckFunctionSignatures
  componentDidUpdate(prevProps: TableProps) {
    const {
      paginationState: { current_page: prevPage }
    } = prevProps;
    const {
      paginationState: { current_page: currentPage },
      onPageChange,
      autoRefresh
    } = this.props;
    const { timer } = this.state;
    if (currentPage && currentPage !== prevPage) onPageChange(currentPage);
    if (!timer && autoRefresh) {
      this.checkForRefresh();
    }
  }

  isRefreshNeeded(): boolean {
    const { autoRefresh, autoRefreshCondition, autoRefreshAction } = this.props;
    if (autoRefresh && autoRefreshCondition && autoRefreshAction) {
      const { paginationState, pagesCollection } = this.props;

      const {
        current_page: currentPage,
        total: totalRows,
        fetching
      } = paginationState;

      const isLoading = fetching || totalRows === null;
      const data =
        !isLoading && has(pagesCollection, currentPage)
          ? pagesCollection[currentPage]
          : [];

      return autoRefreshCondition(data);
    }
    return false;
  }

  checkForRefresh(): void {
    const {
      autoRefresh,
      autoRefreshCondition,
      autoRefreshAction,
      autoRefreshTime
    } = this.props;
    const { timer } = this.state;
    if (timer) {
      clearInterval(timer);
    }
    if (autoRefresh && autoRefreshCondition && autoRefreshAction) {
      let newTimer;
      const doRefresh = () => {
        const {
          paginationState: { current_page: currPage }
        } = this.props;

        if (this.isRefreshNeeded()) {
          autoRefreshAction(currPage || 1);
        } else if (newTimer) {
          console.log('refresh not needed');
          clearInterval(newTimer);
          this.setState({
            timer: undefined
          });
        }
      };
      newTimer = setInterval(doRefresh, autoRefreshTime);
      this.setState({
        timer: newTimer
      });
    }
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
      classes,
      columns,
      paginationState,
      pagesCollection,
      changeSorting,
      requestPage
    } = this.props;

    let { keyField, size, actions, toolbar, title, sortable } = this.props;

    keyField = keyField || 'id';
    size = size || 'small';
    actions = actions || [];
    toolbar = toolbar || [];
    title = title || null;
    sortable = sortable || true;

    const {
      current_page: currentPage,
      per_page: rowsPerPage,
      total: totalRows,
      sorting,
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
        <TableContainer className={classes.container}>
          <TableToolbar
            title={title}
            actions={toolbar}
            state={{ currentPage, rowsPerPage, totalRows, isLoading }}
          />
          {isLoading && (
            <div className={classes.loading}>
              <LinearProgress />
            </div>
          )}
          <Table stickyHeader size={size}>
            <TableHeader
              columns={columns}
              sorting={sorting || {}}
              sortable={sortable}
              changeSorting={changeSorting}
            />
            <TableBody
              data={data}
              keyField={keyField}
              columns={columns}
              actions={actions}
              size={size}
            />
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
