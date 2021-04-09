/* eslint-disable promise/always-return,promise/catch-or-return,react/jsx-props-no-spreading */
// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Box from '@material-ui/core/Box';
import LinearProgress from '@material-ui/core/LinearProgress';
import Icon from '@material-ui/core/Icon';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import PopupState, { bindTrigger, bindMenu } from 'material-ui-popup-state';
import { ConnectTable } from './UI/PaginatedRemoteTable';
import * as JobsActions from '../actions/jobs';
import type { StateType } from '../reducers/types';
import LogsDialog from './UI/LogsDialog';
import * as Api from '../api';
import {
  OUT_TYPE_AR,
  OUT_TYPE_C,
  OUT_TYPE_DESCRIPTION,
  OUT_TYPE_HARMONIZED,
  OUT_TYPE_TRANSCRIPTS
} from '../api/jobs';
import IconButton from './UI/IconButton';
import type { ReadOnlyData, RowActionType } from './UI/Table/types';

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

  openJBrowse = (e, row) => {
    const { pushNotification } = this.props;
    Api.Jobs.openJBrowse(row.id).catch(ex =>
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

  downloadResults = (outputVariable, row) => e => {
    e.preventDefault();
    const { pushNotification } = this.props;
    const jobId = row.id;
    Api.Jobs.genericDownload(
      jobId,
      outputVariable,
      () => {
        const { downloading } = this.state;
        pushNotification('Download started', 'info');
        this.setState({
          downloading: [...downloading, jobId]
        });
      },
      () => {
        const { downloading } = this.state;
        pushNotification('Download completed', 'success');
        this.setState({
          downloading: downloading.filter(i => i !== jobId)
        });
      }
    ).catch(ex => pushNotification(`An error occurred ${ex.message}`, 'error'));
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

  getSaveResultsMenu = (data: ReadOnlyData, size: string) => {
    const { id, status } = data;
    const type = data.output ? data.output.type : null;
    const hasBam = data.output && data.output.outputBamFile;
    const { deletingJobs } = this.props;
    const { downloading } = this.state;
    if (
      deletingJobs.includes(id) ||
      status !== 'completed' ||
      type === OUT_TYPE_C
    ) {
      return null;
    }
    if (downloading.includes(id)) {
      return (
        <IconButton
          key={`action-button-job-${id}`}
          size={size}
          color="inherit"
          onClick={e => e.preventDefault()}
          title="Saving..."
        >
          <Icon className="fas fa-circle-notch fa-spin" fontSize="inherit" />
        </IconButton>
      );
    }
    const items = [['outputFile', 'Save raw output', 'Raw output']];
    if (type === OUT_TYPE_AR)
      items[0] = ['outputFile', 'Save report', 'Report'];
    if (hasBam) {
      items.push(['outputBamFile', 'Save BAM file', 'BAM output']);
    }
    if (OUT_TYPE_HARMONIZED.includes(type))
      items.push([
        'harmonizedFile',
        'Save harmonized output',
        'Harmonized output'
      ]);
    if (OUT_TYPE_TRANSCRIPTS.includes(type))
      items.push([
        'harmonizedTranscriptsFile',
        'Save harmonized transcripts',
        'Harmonized transcripts'
      ]);
    if (OUT_TYPE_DESCRIPTION.includes(type))
      items.push(['description', 'Save metadata file', 'Metadata file']);
    if (items.length === 1) {
      return (
        <IconButton
          key={`action-button-job-${id}`}
          size={size}
          color="inherit"
          onClick={this.downloadResults(items[0][0], data)}
          title={items[0][1]}
        >
          <Icon className="fas fa-save" fontSize="inherit" />
        </IconButton>
      );
    }
    return (
      <PopupState
        variant="popover"
        popupId={`popup-menu-job-${id}`}
        key={`popup-menu-job-${id}`}
      >
        {popupState => {
          const iconProps = {
            size,
            color: 'inherit',
            title: 'Save',
            // $FlowFixMe: flow is stupid
            ...bindTrigger(popupState),
            key: `popup-button-job-${id}`
          };
          return (
            <>
              <IconButton {...iconProps}>
                <Icon className="fas fa-save" fontSize="inherit" />
              </IconButton>
              <Menu {...bindMenu(popupState)} key={`popup-menu-job-${id}`}>
                {items.map(i => (
                  <MenuItem
                    key={`popup-menu-item-${i[0]}-job-${id}`}
                    onClick={e => {
                      popupState.close();
                      this.downloadResults(i[0], data)(e);
                    }}
                  >
                    {i[2]}
                  </MenuItem>
                ))}
              </Menu>
            </>
          );
        }}
      </PopupState>
    );
  };

  getActions(): RowActionType[] {
    const { deletingJobs, submittingJobs } = this.props;
    const { downloading } = this.state;

    const cd = r => deletingJobs.includes(r.id);
    const cs = r => submittingJobs.includes(r.id);
    const cw = r => downloading.includes(r.id);
    const isReady = r => r.status === 'ready';
    const isCompleted = r => r.status === 'completed';
    const isConfirmation = r => r.output && r.output.type === OUT_TYPE_C;
    const isReport = r => r.output && r.output.type === OUT_TYPE_AR;
    const hasBrowser = r =>
      r.output &&
      r.output.outputJBrowseFile &&
      typeof r.output.outputJBrowseFile === 'object';

    return [
      {
        shown: r => !cd(r) && isReady(r) && !cs(r),
        icon: 'fas fa-play',
        color: 'primary',
        onClick: this.handleSubmitJob,
        tooltip: 'Submit'
      },
      {
        disabled: true,
        shown: r => !cd(r) && isReady(r) && cs(r),
        icon: 'fas fa-circle-notch fa-spin',
        color: 'primary',
        tooltip: 'Submitting...'
      },
      {
        shown: r =>
          !cd(r) && ['processing', 'completed', 'failed'].includes(r.status),
        icon: 'fas fa-file-alt',
        tooltip: 'Logs',
        onClick: this.handleLogsSelectJob
      },
      this.getSaveResultsMenu,
      {
        shown: r =>
          !cd(r) &&
          isCompleted(r) &&
          Api.Settings.isLocal() &&
          !isConfirmation(r),
        icon: 'fas fa-folder-open',
        tooltip: 'Open results folder',
        onClick: this.openResultsFolder
      },
      {
        shown: r => !cd(r) && isCompleted(r) && !cw(r) && hasBrowser(r),
        icon: 'fas fa-dna',
        tooltip: 'Show Genome Browser',
        onClick: this.openJBrowse
      },
      {
        shown: r => !cd(r) && isCompleted(r) && !cw(r) && isReport(r),
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
