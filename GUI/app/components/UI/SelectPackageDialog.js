// @flow
/* eslint-disable no-nested-ternary */
import * as React from 'react';
import Dialog from '@material-ui/core/Dialog';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import useMediaQuery from '@material-ui/core/useMediaQuery';
import { useTheme } from '@material-ui/core/styles';
import DialogActions from '@material-ui/core/DialogActions';
import Button from '@material-ui/core/Button';
import { Form, Formik } from 'formik';
import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import * as Api from '../../api';
import TableField from '../Form/TableField';
import { pushNotificationSimple } from '../../actions/notifications';
import type { PushNotificationFunction } from '../../types/notifications';

type Props = {
  open: boolean,
  onClose: () => void,
  pushNotification: PushNotificationFunction
};

function SelectPackageDialog({ open, onClose, pushNotification }: Props) {
  if (!Api.Settings.isLocal() || !Api.Settings.isConfigured()) return null;
  const theme = useTheme();
  const [selectedPackage, setSelectedPackage] = React.useState<?string>(null);
  const fullScreen = useMediaQuery(theme.breakpoints.down('md'));
  const clickInstallPackage = () => {
    if (selectedPackage) {
      console.log(selectedPackage);
    }
  };
  return (
    <Dialog fullScreen={fullScreen} open={open} onClose={onClose}>
      <DialogTitle>Install packages</DialogTitle>
      <DialogContent>
        <Formik
          initialValues={{
            package: ''
          }}
          onSubmit={() => undefined}
        >
          <TableField
            name="source_sample_group"
            required
            single
            getData={() => Api.Docker.listPackages()}
            onError={e =>
              pushNotification(`An error occurred: ${e.message}!`, 'error')
            }
            onChange={v => setSelectedPackage(v)}
            keyField="name"
            label="Select one package"
            columns={[
              {
                dataField: 'title',
                label: 'Title'
              },
              {
                dataField: 'description',
                label: 'Description'
              }
            ]}
          />
        </Formik>
      </DialogContent>
      <DialogActions>
        <Button onClick={clickInstallPackage} color="primary" autoFocus>
          Install
        </Button>
      </DialogActions>
    </Dialog>
  );
}

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      pushNotification: pushNotificationSimple
    },
    dispatch
  );
}

// $FlowFixMe: Flow is stupid
export default connect(() => ({}), mapDispatchToProps)(SelectPackageDialog);
