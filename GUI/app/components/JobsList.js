/* eslint-disable promise/always-return,promise/catch-or-return */
// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Box from '@material-ui/core/Box';
import LinearProgress from '@material-ui/core/LinearProgress';
import { ConnectTable } from './UI/PaginatedRemoteTable';
import * as JobsActions from '../actions/jobs';
import type { StateType } from '../reducers/types';
import LogsDialog from './UI/LogsDialog';
import * as Api from '../api';

const JobsTable = ConnectTable(
  (state: StateType) => ({
    paginationState: state.jobs.jobsList.state,
    pagesCollection: state.jobs.jobsList.pages
  }),
  {
    changeRowsPerPage: JobsActions.setPerPage,
    changeSorting: JobsActions.setSorting,
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
  deletingJobs: number[],
  submittingJobs: number[],
  submitJob: (number, ?number) => void,
  refreshPage: number => void,
  deleteJob: number => void,
  pushNotification: (string, 'success' | 'warning' | 'error' | 'info') => void,
  classes: {
    root: *,
    loading: *
  }
};

type State = {
  isLoading: boolean,
  logsOpen: boolean,
  logsSelectedJobId: ?number,
  currentPage: ?number,
  downloading: number[]
};

class JobsList extends React.Component<Props, State> {
  constructor(props) {
    super(props);
    this.state = {
      isLoading: false,
      logsOpen: false,
      logsSelectedJobId: null,
      currentPage: null,
      downloading: []
    };
  }

  handleLogsClose = () => {
    const { currentPage } = this.state;
    const { refreshPage } = this.props;
    this.setState({ logsOpen: false, logsSelectedJobId: null });
    if (currentPage) refreshPage(currentPage);
  };

  handleJobDelete = (e, row) => {
    const { deleteJob } = this.props;
    deleteJob(row.id);
    e.preventDefault();
  };

  openResultsFolder = (e, row) => {
    const { pushNotification } = this.props;
    Api.Jobs.openLocalFolder(row.id).catch(ex =>
      pushNotification(`An error occurred ${ex.message}`, 'error')
    );
    e.preventDefault();
  };

  openReport = (e, row) => {
    const { pushNotification } = this.props;
    Api.Jobs.openReport(row.id).catch(ex =>
      pushNotification(`An error occurred ${ex.message}`, 'error')
    );
    e.preventDefault();
  };

  downloadResults = (e, row) => {
    const { pushNotification } = this.props;
    const jobId = row.id;
    Api.Jobs.download(
      jobId,
      () => {
        const { downloading } = this.state;
        this.setState({
          downloading: [...downloading, jobId]
        });
      },
      () => {
        const { downloading } = this.state;
        this.setState({
          downloading: downloading.filter(i => i !== jobId)
        });
      }
    ).catch(ex => pushNotification(`An error occurred ${ex.message}`, 'error'));
    e.preventDefault();
  };

  handleLogsSelectJob = (e, row) => {
    this.setState({
      logsOpen: true,
      logsSelectedJobId: row.id
    });
    e.preventDefault();
  };

  handlePageChange = (currentPage: number) => {
    this.setState({ currentPage });
  };

  handleSubmitJob = (e, row) => {
    const { submitJob } = this.props;
    const { currentPage } = this.state;
    submitJob(row.id, currentPage);
    e.preventDefault();
  };

  handleRefresh = (e, state) => {
    if (!state.isLoading) {
      const { refreshPage } = this.props;
      refreshPage(state.currentPage || 1);
    }
    e.preventDefault();
  };

  getActions() {
    const { deletingJobs, submittingJobs } = this.props;
    const { downloading } = this.state;

    const cd = r => deletingJobs.includes(r.id);
    const cs = r => submittingJobs.includes(r.id);
    const cw = r => downloading.includes(r.id);
    const isReport = r => r.type === 'diff_expr_analysis_job_type';

    return [
      {
        shown: r => !cd(r) && r.status === 'ready' && !cs(r),
        icon: 'fas fa-play',
        color: 'primary',
        onClick: this.handleSubmitJob,
        tooltip: 'Submit'
      },
      {
        disabled: true,
        shown: r => !cd(r) && r.status === 'ready' && cs(r),
        icon: 'fas fa-circle-notch fa-spin',
        color: 'primary',
        tooltip: 'Submitting...'
      },
      {
        shown: r => !cd(r) && r.status !== 'ready' && r.status !== 'queued',
        icon: 'fas fa-file-alt',
        tooltip: 'Logs',
        onClick: this.handleLogsSelectJob
      },
      {
        shown: r => !cd(r) && r.status === 'completed' && cw(r),
        icon: 'fas fa-circle-notch fa-spin',
        tooltip: 'Saving...'
      },
      {
        shown: r => !cd(r) && r.status === 'completed' && !cw(r),
        icon: 'fas fa-save',
        tooltip: 'Save results',
        onClick: this.downloadResults
      },
      {
        shown: r =>
          !cd(r) && r.status === 'completed' && Api.Settings.isLocal(),
        icon: 'fas fa-folder-open',
        tooltip: 'Open results folder',
        onClick: this.openResultsFolder
      },
      {
        shown: r => !cd(r) && r.status === 'completed' && !cw(r) && isReport(r),
        icon: 'fas fa-eye',
        tooltip: 'Show report',
        onClick: this.openReport
      },
      {
        shown: r => !cd(r) && r.status !== 'processing',
        color: 'secondary',
        icon: 'fas fa-trash',
        tooltip: 'Delete',
        onClick: this.handleJobDelete
      }
    ];
  }

  getToolbarActions() {
    return [
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
    const { deletingJobs, classes } = this.props;
    const { isLoading, logsOpen, logsSelectedJobId } = this.state;

    return (
      <>
        <Box className={classes.root}>
          {isLoading && (
            <div className={classes.loading}>
              <LinearProgress />
            </div>
          )}
          <div>
            <JobsTable
              title="Jobs list"
              onPageChange={this.handlePageChange}
              toolbar={this.getToolbarActions()}
              actions={this.getActions()}
              columns={[
                {
                  dataField: 'name',
                  label: 'Name'
                },
                {
                  dataField: 'readable_type',
                  sortingField: 'job_type',
                  label: 'Type'
                },
                {
                  dataField: 'status',
                  label: 'Status',
                  format: (v, r) => {
                    if (deletingJobs.includes(r.id)) return 'Deleting';
                    return Api.Utils.capitalize(v);
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
