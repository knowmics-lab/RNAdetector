import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import JobsList from '../components/JobsList';
import * as JobsActions from '../actions/jobs';

function mapStateToProps() {
  return {};
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
