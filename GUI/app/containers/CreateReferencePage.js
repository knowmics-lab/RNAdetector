import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import { push } from 'connected-react-router';
import * as ReferencesActions from '../actions/references';
import CreateReference from '../components/CreateReference';

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      push
    },
    dispatch
  );
}

export default connect(() => ({}), mapDispatchToProps)(CreateReference);
