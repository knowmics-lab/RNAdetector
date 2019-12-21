import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import * as JobsActions from '../actions/jobs';
import JobsList from '../components/JobsList';

function mapStateToProps(state) {
  return {
    ...state.jobs.jobsList
  };
}

function mapDispatchToProps(dispatch) {
  return bindActionCreators(JobsActions, dispatch);
}

export default connect(mapStateToProps, mapDispatchToProps)(JobsList);
