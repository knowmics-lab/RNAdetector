import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import JobsList from '../components/JobsList';
import * as JobsActions from '../actions/jobs';

function mapStateToProps(state) {
  return {
    submittingJobs: state.jobs.jobs.submitting,
    deletingJobs: state.jobs.jobs.deleting
  };
}

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      submitJob: JobsActions.submitJob,
      refreshPage: JobsActions.refreshPage,
      deleteJob: JobsActions.deleteJob
    },
    dispatch
  );
}

export default connect(mapStateToProps, mapDispatchToProps)(JobsList);
