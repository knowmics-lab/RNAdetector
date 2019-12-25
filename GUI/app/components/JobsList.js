/* eslint-disable promise/always-return,promise/catch-or-return */
// @flow
import React, { Component } from 'react';
import { withStyles } from '@material-ui/core/styles';
import { Typography, Paper, Box, Icon } from '@material-ui/core';
import LinearProgress from '@material-ui/core/LinearProgress';
import { ConnectTable } from './UI/PaginatedRemoteTable';
import * as JobsActions from '../actions/jobs';
import type { StateType } from '../reducers/types';
import IconButton from './UI/IconButton';
import LogsDialog from './UI/LogsDialog';

const JobsTable = ConnectTable(
  (state: StateType) => ({
    paginationState: state.jobs.jobsList.state,
    pagesCollection: state.jobs.jobsList.pages
  }),
  {
    changeRowsPerPage: JobsActions.setPerPage,
    resetErrorMessage: JobsActions.jobsListResetErrorMessage,
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
  isLoading: boolean,
  logsOpen: boolean,
  logsSelectedJobId: ?number
};

class JobsList extends Component<Props, State> {
  constructor(props) {
    super(props);
    this.state = {
      isLoading: false,
      logsOpen: false,
      logsSelectedJobId: null
    };
  }

  handleLogsClose = () =>
    this.setState({ logsOpen: false, logsSelectedJobId: null });

  handleLogsSelectJob = (jobId: number) =>
    this.setState({
      logsOpen: true,
      logsSelectedJobId: jobId
    });

  render() {
    const { classes } = this.props;
    const { isLoading, logsOpen, logsSelectedJobId } = this.state;
    return (
      <>
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
                    id: 'name',
                    label: 'Name'
                  },
                  {
                    id: 'readable_type',
                    label: 'Type'
                  },
                  {
                    id: 'status',
                    label: 'Status'
                  },
                  {
                    id: 'created_at_diff',
                    label: 'Created at'
                  },
                  {
                    id: 'id',
                    align: 'center',
                    label: 'Action',
                    format: row => {
                      // eslint-disable-next-line radix
                      // 'ready' | 'queued' | 'processing' | 'completed' | 'failed'
                      const components = [];
                      if (row.status === 'ready') {
                        components.push(
                          <IconButton
                            title="Submit"
                            color="primary"
                            href="/href"
                            key={`${row.id}-submit`}
                          >
                            <Icon className="fas fa-play" />
                          </IconButton>
                        );
                      }
                      if (row.status !== 'ready' && row.status !== 'queued') {
                        components.push(
                          <IconButton
                            title="Logs"
                            onClick={() => this.handleLogsSelectJob(row.id)}
                            key={`${row.id}-logs`}
                          >
                            <Icon className="fas fa-file-alt" />
                          </IconButton>
                        );
                      }
                      if (row.status === 'completed') {
                        components.push(
                          <IconButton
                            title="Save results"
                            href="/test"
                            key={`${row.id}-save`}
                          >
                            <Icon className="fas fa-save" />
                          </IconButton>
                        );
                      }
                      if (row.status !== 'processing') {
                        components.push(
                          <IconButton
                            title="Delete"
                            color="secondary"
                            href="/test"
                            key={`${row.id}-delete`}
                          >
                            <Icon className="fas fa-trash" />
                          </IconButton>
                        );
                      }
                      return <>{components}</>;
                    }
                  }
                ]}
              />
            </div>
          </Paper>
        </Box>
        <LogsDialog
          jobId={logsSelectedJobId}
          open={logsOpen}
          onClose={this.handleLogsClose}
        />
      </>
    );
  }
}

export default withStyles(style)(JobsList);
