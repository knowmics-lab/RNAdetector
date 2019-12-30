import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import ReferencesList from '../components/ReferencesList';
import * as ReferencesActions from '../actions/references';

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      refreshPage: ReferencesActions.refreshPage,
      deleteReference: ReferencesActions.deleteReference
    },
    dispatch
  );
}

export default connect(() => ({}), mapDispatchToProps)(ReferencesList);
