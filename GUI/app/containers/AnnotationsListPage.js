import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import AnnotationsList from '../components/AnnotationsList';
import * as AnnotationsActions from '../actions/annotations';

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      refreshPage: AnnotationsActions.refreshPage,
      deleteAnnotation: AnnotationsActions.deleteAnnotation
    },
    dispatch
  );
}

export default connect(() => ({}), mapDispatchToProps)(AnnotationsList);
