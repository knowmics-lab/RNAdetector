// @flow
import React, { Component } from 'react';
import { withStyles } from '@material-ui/core/styles';
import { Typography, Paper, Box } from '@material-ui/core';
import moment from 'moment';
import { ConnectTable } from './UI/PaginatedRemoteTable';
import * as JobsActions from '../actions/jobs';
import type { StateType } from '../reducers/types';

const JobsTable = ConnectTable(
  (state: StateType) => ({
    paginationState: state.jobs.jobsList.state,
    pagesCollection: state.jobs.jobsList.pages
  }),
  {
    changeRowsPerPage: JobsActions.setPerPage,
    handleSnackbarClose: JobsActions.jobsListResetLoading,
    requestPage: JobsActions.requestPage
  }
);

const style = theme => ({
  root: {
    padding: theme.spacing(3, 2)
  }
});

type Props = {
  classes: {
    root: *
  }
};

class JobsList extends Component<Props> {
  props: Props;

  render() {
    const { classes } = this.props;
    return (
      <Box>
        <Paper className={classes.root}>
          <Typography variant="h5" component="h3">
            Jobs list
          </Typography>
          <Typography component="p" />
          <div>
            <JobsTable
              columns={[
                {
                  id: 'type',
                  label: 'Type',
                  minWidth: 100
                },
                {
                  id: 'status',
                  label: 'Status',
                  minWidth: 100
                },
                {
                  id: 'created_at',
                  label: 'Created at',
                  format: value => {
                    return moment(value).format('YYYY-MM-DD HH:mm:ss');
                  },
                  minWidth: 150
                }
              ]}
            />
          </div>
        </Paper>
      </Box>
    );
  }
}

export default withStyles(style)(JobsList);
