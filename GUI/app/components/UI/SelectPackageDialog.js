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
import { Formik } from 'formik';
import { bindActionCreators } from 'redux';
import { connect } from 'react-redux';
import { Tooltip } from '@material-ui/core';
import { push as redirectAction } from 'connected-react-router';
import * as Api from '../../api';
import TableField from '../Form/TableField';
import { pushNotificationSimple } from '../../actions/notifications';
import type { PushNotificationFunction } from '../../types/notifications';
import { refreshPage as refreshJobsAction } from '../../actions/jobs';
import { JOBS } from '../../constants/routes.json';

type Props = {
  open: boolean,
  onClose: () => void,
  pushNotification: PushNotificationFunction,
  refreshJobs: () => void,
  redirect: mixed => void
};

function SelectPackageDialog({
  open,
  onClose,
  pushNotification,
  refreshJobs,
  redirect
}: Props) {
  if (!Api.Settings.isConfigured()) return null;
  const theme = useTheme();
  const [isInstalling, setIsInstalling] = React.useState<boolean>(false);
  const [selectedPackage, setSelectedPackage] = React.useState<?(string[])>(
    null
  );

  const fullScreen = useMediaQuery(theme.breakpoints.down('md'));
  const clickInstallPackage = async () => {
    if (selectedPackage && selectedPackage.length > 0) {
      setIsInstalling(true);
      const data = await Api.References.install(selectedPackage);
      const { data: job } = data;
      await Api.Jobs.submitJob(job.id);
      refreshJobs();
      pushNotification('An installation request has been created', 'success');
      redirect(JOBS);
    } else {
      pushNotification(
        'You must select at least one package to install',
        'error'
      );
    }
  };
  return (
    <>
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
              name="package"
              required
              getData={() => Api.References.listPackages()}
              onError={e =>
                pushNotification(`An error occurred: ${e.message}!`, 'error')
              }
              onChange={v => setSelectedPackage(v)}
              keyField="name"
              label="Select packages"
              columns={[
                {
                  dataField: 'title',
                  label: 'Title'
                },
                {
                  dataField: 'description',
                  label: 'Description'
                },
                {
                  dataField: 'needsUpdate',
                  label: '',
                  format: v => {
                    return v ? (
                      <Tooltip title="An update is available!">
                        <i className="fas fa-exclamation-circle" />
                      </Tooltip>
                    ) : (
                      ''
                    );
                  }
                }
              ]}
            />
          </Formik>
        </DialogContent>
        <DialogActions>
          <Button
            color="secondary"
            onClick={clickInstallPackage}
            disabled={isInstalling}
            autoFocus
          >
            {isInstalling ? 'Please wait...' : 'Install'}
          </Button>
          <Button onClick={onClose} color="primary" autoFocus>
            Close
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
}

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      pushNotification: pushNotificationSimple,
      refreshJobs: refreshJobsAction,
      redirect: redirectAction
    },
    dispatch
  );
}

// $FlowFixMe: Flow is stupid
export default connect(() => ({}), mapDispatchToProps)(SelectPackageDialog);
