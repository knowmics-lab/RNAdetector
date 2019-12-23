/* eslint-disable promise/always-return,promise/catch-or-return */
// @flow
import React, { Component } from 'react';
import { withStyles } from '@material-ui/core/styles';
import { Typography, Paper, Box, Icon } from '@material-ui/core';
import moment from 'moment';
import LinearProgress from '@material-ui/core/LinearProgress';
import { ConnectTable } from './UI/PaginatedRemoteTable';
import * as JobsActions from '../actions/jobs';
import * as Api from '../api';
import type { StateType } from '../reducers/types';
import Button from './UI/Button';
import IconButton from "./UI/IconButton";

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
  },
  loading: {
    width: '100%',
    '& > * + *': {
      marginTop: theme.spacing(2)
    }
  }
});

type Props = {
  classes: {
    root: *,
    loading: *
  }
};

type State = {
  isLoading: boolean
};

class JobsList extends Component<Props, State> {
  constructor(props) {
    super(props);
    this.state = {
      isLoading: false
    };
  }

  render() {
    const { classes } = this.props;
    const { isLoading } = this.state;
    return (
      <Box>
        <Paper className={classes.root}>
          <Typography variant="h5" component="h3">
            Jobs list
          </Typography>
          <Typography component="p" />
          {isLoading && (
            <div className={classes.loading}>
              <LinearProgress />
            </div>
          )}
          <div>
            <JobsTable
              columns={[
                {
                  id: 'type',
                  label: 'Type',
                  minWidth: 100,
                  format: value => {
                    return Api.Utils.capitalize(
                      Api.Utils.dashToWordString(value)
                        .replace(' Job Type', '')
                        .replace('Rna', 'RNA')
                        .replace('Dna', 'DNA')
                    );
                  }
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
                },
                {
                  id: 'id',
                  label: 'Action',
                  format: value => {
                    // eslint-disable-next-line radix
                    return (
                      <>
                        <IconButton title="Logs" href="/test">
                          <Icon className="fas fa-file-alt" />
                        </IconButton>
                        <IconButton
                          title="Delete"
                          color="secondary"
                          href="/test"
                        >
                          <Icon className="fas fa-trash" />
                        </IconButton>
                      </>
                    );
                  }
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
