/* eslint-disable promise/always-return,promise/catch-or-return */
// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Box from '@material-ui/core/Box';
import capitalize from '@material-ui/core/utils/capitalize';
import { ConnectTable } from './UI/PaginatedRemoteTable';
import * as ReferencesActions from '../actions/references';
import type { StateType } from '../reducers/types';

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
  classes: {
    root: *,
    loading: *
  }
};

type State = {
  currentPage: ?number
};

class ReferencesList extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.state = {
      currentPage: null
    };
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
    return [
      {
        align: 'right',
        shown: true,
        icon: 'fas fa-plus',
        tooltip: 'Add',
        onClick: () => console.log('TODO')
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

  render() {
    const { classes } = this.props;

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
                    return Object.keys(value)
                      .map(capitalize)
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
      </Box>
    );
  }
}

export default withStyles(style)(ReferencesList);
