import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import { push } from 'connected-react-router';
import { refreshPage } from '../actions/jobs';
import { pushNotificationSimple } from '../actions/notifications';
import CreateAnnotation from '../components/CreateAnnotation';

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      pushNotification: pushNotificationSimple,
      redirect: push,
      refreshJobs: refreshPage
    },
    dispatch
  );
}

export default connect(() => ({}), mapDispatchToProps)(CreateAnnotation);
