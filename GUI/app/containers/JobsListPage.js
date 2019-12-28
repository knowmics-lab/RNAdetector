import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import JobsList from '../components/JobsList';
import * as JobsActions from '../actions/jobs';

function mapStateToProps(state) {
  return {
    submittingJobs: state.jobs.jobs.submitting
  };
}

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      submitJob: JobsActions.submitJob,
      refreshPage: JobsActions.refreshPage
    },
    dispatch
  );
}

export default connect(mapStateToProps, mapDispatchToProps)(JobsList);
