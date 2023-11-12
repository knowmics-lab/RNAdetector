import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import { push } from 'connected-react-router';
import Home from '../components/Home';

function mapStateToProps(state) {
  return {};
}

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      redirect: push
    },
    dispatch
  );
}

export default connect(mapStateToProps, mapDispatchToProps)(Home);
