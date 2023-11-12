import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import { push } from 'connected-react-router';
import ReferencesList from '../components/ReferencesList';
import * as ReferencesActions from '../actions/references';

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      refreshPage: ReferencesActions.refreshPage,
      deleteReference: ReferencesActions.deleteReference,
      push
    },
    dispatch
  );
}

export default connect(() => ({}), mapDispatchToProps)(ReferencesList);
