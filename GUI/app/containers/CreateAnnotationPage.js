import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import { push } from 'connected-react-router';
import { refreshPage as refreshJobs } from '../actions/jobs';
import { refreshPage as refreshAnnotations } from '../actions/annotations';
import { pushNotificationSimple } from '../actions/notifications';
import CreateAnnotation from '../components/CreateAnnotation';

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      pushNotification: pushNotificationSimple,
      redirect: push,
      refreshJobs,
      refreshAnnotations
    },
    dispatch
  );
}

export default connect(() => ({}), mapDispatchToProps)(CreateAnnotation);
