/* eslint-disable promise/always-return,promise/catch-or-return */
// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Box from '@material-ui/core/Box';
import capitalize from '@material-ui/core/utils/capitalize';
import { ConnectTable } from './UI/PaginatedRemoteTable';
import * as ReferencesActions from '../actions/references';
import * as Api from '../api';
import type { StateType } from '../reducers/types';
import { CREATE_REFERENCE } from '../constants/routes';
import SelectPackageDialog from './UI/SelectPackageDialog';

const ReferencesTable = ConnectTable(
  (state: StateType) => ({
    paginationState: state.references.referencesList.state,
    pagesCollection: state.references.referencesList.pages
  }),
  {
    changeRowsPerPage: ReferencesActions.setPerPage,
    changeSorting: ReferencesActions.setSorting,
    requestPage: ReferencesActions.requestPage
  }
);

const style = theme => ({
  root: {
    padding: theme.spacing(3, 2)
  },
  loading: {
    width: '100%',
    '& > * + *': {
      marginTop: theme.spacing(2)
    }
  }
});

type Props = {
  refreshPage: number => void,
  deleteReference: (number, ?number) => void,
  push: (*) => void,
  classes: {
    root: *,
    loading: *
  }
};

type State = {
  currentPage: ?number,
  selectDialogOpen: boolean,
  isOnline: boolean
};

class ReferencesList extends React.Component<Props, State> {
  props: Props;

  timer: *;

  constructor(props) {
    super(props);
    this.state = {
      currentPage: null,
      selectDialogOpen: false,
      isOnline: false
    };
  }

  componentDidMount(): void {
    this.timer = setInterval(
      () =>
        Api.Utils.isOnline().then(r => {
          this.setState({
            isOnline: r
          });
        }),
      5000
    );
  }

  componentWillUnmount(): void {
    if (this.timer) {
      clearInterval(this.timer);
    }
  }

  handleReferenceDelete = (e, row) => {
    const { currentPage } = this.state;
    const { deleteReference } = this.props;
    deleteReference(row.id, currentPage);
    e.preventDefault();
  };

  handlePageChange = (currentPage: number) => this.setState({ currentPage });

  handleRefresh = (e, state) => {
    if (!state.isLoading) {
      const { refreshPage } = this.props;
      refreshPage(state.currentPage || 1);
    }
    e.preventDefault();
  };

  getActions() {
    return [
      {
        shown: true,
        color: 'secondary',
        icon: 'fas fa-trash',
        tooltip: 'Delete',
        onClick: this.handleReferenceDelete
      }
    ];
  }

  getToolbarActions() {
    const { isOnline } = this.state;
    return [
      {
        align: 'right',
        shown: true,
        icon: 'fas fa-plus',
        tooltip: 'Add',
        onClick: () => {
          const { push } = this.props;
          push(CREATE_REFERENCE);
        }
      },
      {
        align: 'right',
        shown:
          Api.Settings.isLocal() && Api.Settings.isConfigured() && isOnline,
        icon: 'fas fa-download',
        tooltip: 'Install from repository',
        onClick: () =>
          this.setState({
            selectDialogOpen: true
          })
      },
      {
        align: 'right',
        shown: true,
        icon: 'fas fa-redo',
        disabled: s => s.isLoading,
        tooltip: 'Refresh',
        onClick: this.handleRefresh
      }
    ];
  }

  closeSelectDialog = () => {
    this.setState({
      selectDialogOpen: false
    });
  };

  render() {
    const { classes } = this.props;
    const { selectDialogOpen } = this.state;

    return (
      <Box className={classes.root}>
        <div>
          <ReferencesTable
            title="References list"
            onPageChange={this.handlePageChange}
            toolbar={this.getToolbarActions()}
            actions={this.getActions()}
            columns={[
              {
                dataField: 'name',
                label: 'Name'
              },
              {
                dataField: 'available_for',
                label: 'Indexed for',
                disableSorting: true,
                format: value => {
                  if (typeof value === 'object' && value) {
                    return Object.entries(value)
                      .filter(([, v]) => v)
                      .map(([k]) => capitalize(k))
                      .join(', ');
                  }
                  return '';
                }
              },
              {
                dataField: 'created_at_diff',
                sortingField: 'created_at',
                label: 'Created at'
              },
              'actions'
            ]}
          />
        </div>
        <SelectPackageDialog
          open={selectDialogOpen}
          onClose={this.closeSelectDialog}
        />
      </Box>
    );
  }
}

export default withStyles(style)(ReferencesList);
