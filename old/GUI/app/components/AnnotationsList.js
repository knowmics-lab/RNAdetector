/* eslint-disable promise/always-return,promise/catch-or-return */
// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Box from '@material-ui/core/Box';
import { ConnectTable } from './UI/PaginatedRemoteTable';
import * as AnnotationsActions from '../actions/annotations';
import type { StateType } from '../reducers/types';
import { CREATE_ANNOTATION } from '../constants/routes.json';

const AnnotationsTable = ConnectTable(
  (state: StateType) => ({
    paginationState: state.annotations.annotationsList.state,
    pagesCollection: state.annotations.annotationsList.pages
  }),
  {
    changeRowsPerPage: AnnotationsActions.setPerPage,
    changeSorting: AnnotationsActions.setSorting,
    requestPage: AnnotationsActions.requestPage
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
  deleteAnnotation: (number, ?number) => void,
  push: (*) => void,
  classes: {
    root: *,
    loading: *
  }
};

type State = {
  currentPage: ?number
};

class AnnotationsList extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.state = {
      currentPage: null
    };
  }

  handleAnnotationDelete = (e, row) => {
    const { currentPage } = this.state;
    const { deleteAnnotation } = this.props;
    deleteAnnotation(row.id, currentPage);
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
        onClick: this.handleAnnotationDelete
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
        onClick: () => {
          const { push } = this.props;
          push(CREATE_ANNOTATION);
        }
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
          <AnnotationsTable
            title="Annotations list"
            onPageChange={this.handlePageChange}
            toolbar={this.getToolbarActions()}
            actions={this.getActions()}
            columns={[
              {
                dataField: 'name',
                label: 'Name'
              },
              {
                dataField: 'type',
                label: 'Type',
                format: value =>
                  typeof value === 'string' ? value.toUpperCase() : ''
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

export default withStyles(style)(AnnotationsList);
