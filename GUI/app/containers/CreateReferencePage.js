import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import { push } from 'connected-react-router';
import { refreshPage as refreshJobs } from '../actions/jobs';
import { refreshPage as refreshReferences } from '../actions/references';
import CreateReference from '../components/CreateReference';
import { pushNotificationSimple } from '../actions/notifications';

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      pushNotification: pushNotificationSimple,
      redirect: push,
      refreshJobs,
      refreshReferences
    },
    dispatch
  );
}

export default connect(() => ({}), mapDispatchToProps)(CreateReference);
