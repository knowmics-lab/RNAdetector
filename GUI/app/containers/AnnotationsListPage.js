import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import { push } from 'connected-react-router';
import AnnotationsList from '../components/AnnotationsList';
import * as AnnotationsActions from '../actions/annotations';

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      refreshPage: AnnotationsActions.refreshPage,
      deleteAnnotation: AnnotationsActions.deleteAnnotation,
      push
    },
    dispatch
  );
}

export default connect(() => ({}), mapDispatchToProps)(AnnotationsList);
