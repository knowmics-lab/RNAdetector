import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import { push } from 'connected-react-router';
import { refreshPage as refreshJobs } from '../../actions/jobs';
import { pushNotificationSimple } from '../../actions/notifications';
import DiffExpr from '../../components/Analysis/DiffExpr';

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      pushNotification: pushNotificationSimple,
      redirect: push,
      refreshJobs
    },
    dispatch
  );
}

export default connect(() => ({}), mapDispatchToProps)(DiffExpr);
