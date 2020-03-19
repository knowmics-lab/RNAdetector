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
import * as Api from '../../api';
import TableField from '../Form/TableField';
import { pushNotificationSimple } from '../../actions/notifications';
import type { PushNotificationFunction } from '../../types/notifications';
import { refreshPage as refreshAnnotationsAction } from '../../actions/annotations';
import * as ReferencesActions from '../../actions/references';

type InstallingProps = {
  open: boolean,
  content: string,
  closeButton: boolean,
  requestClose: () => void
};

function InstallingDialog({
  open,
  content,
  closeButton,
  requestClose
}: InstallingProps) {
  const logRef = React.createRef();
  React.useEffect(() => {
    if (logRef.current) {
      logRef.current.scrollIntoView({ behavior: 'smooth' });
    }
  });

  return (
    <Dialog fullScreen open={open}>
      <DialogTitle>Installing package</DialogTitle>
      <DialogContent>
        <pre>{content}</pre>
        <div ref={logRef} />
      </DialogContent>
      {closeButton && (
        <DialogActions>
          <Button onClick={requestClose} color="primary" autoFocus>
            Close
          </Button>
        </DialogActions>
      )}
    </Dialog>
  );
}

type Props = {
  open: boolean,
  onClose: () => void,
  pushNotification: PushNotificationFunction,
  refreshAnnotations: () => void,
  refreshPage: number => void
};

type LogState = {
  open: boolean,
  content: string,
  closeButton: boolean
};

function SelectPackageDialog({
  open,
  onClose,
  pushNotification,
  refreshAnnotations,
  refreshPage
}: Props) {
  if (!Api.Settings.isLocal() || !Api.Settings.isConfigured()) return null;
  const theme = useTheme();
  const [selectedPackage, setSelectedPackage] = React.useState<?string>(null);
  const [logState, setLogState] = React.useState<LogState>({
    open: false,
    content: '',
    closeButton: false
  });
  const onCloseRequest = () =>
    setLogState(prev => ({
      ...prev,
      open: false
    }));

  const fullScreen = useMediaQuery(theme.breakpoints.down('md'));
  const clickInstallPackage = () => {
    if (selectedPackage) {
      setLogState({
        open: true,
        content: '',
        closeButton: false
      });
      Api.Docker.installPackage(
        selectedPackage,
        m =>
          setLogState(prev => ({
            ...prev,
            content: `${prev.content.trimEnd()}\n${m}`
          })),
        e =>
          setLogState(prev => ({
            ...prev,
            content: `${prev.content.trimEnd()}\nAn error occurred: ${
              e.message
            }\n`,
            closeButton: true
          })),
        c => {
          const newContent =
            c === 0 ? 'Package installed!' : 'An unknown error occurred!';
          refreshPage(1);
          refreshAnnotations();
          setTimeout(
            () =>
              setLogState(prev => ({
                ...prev,
                content: `${prev.content.trimEnd()}\n${newContent}\n`,
                closeButton: true
              })),
            3000
          );
        }
      );
    } else {
      pushNotification('You must select one package to install', 'error');
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
          <Button onClick={clickInstallPackage} color="secondary" autoFocus>
            Install
          </Button>
          <Button onClick={onClose} color="primary" autoFocus>
            Close
          </Button>
        </DialogActions>
      </Dialog>
      <InstallingDialog
        open={logState.open}
        content={logState.content}
        closeButton={logState.closeButton}
        requestClose={onCloseRequest}
      />
    </>
  );
}

function mapDispatchToProps(dispatch) {
  return bindActionCreators(
    {
      pushNotification: pushNotificationSimple,
      refreshAnnotations: refreshAnnotationsAction,
      refreshPage: ReferencesActions.refreshPage
    },
    dispatch
  );
}

// $FlowFixMe: Flow is stupid
export default connect(() => ({}), mapDispatchToProps)(SelectPackageDialog);
