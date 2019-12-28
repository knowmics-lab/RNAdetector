// @flow
import * as React from 'react';
import { connect } from 'react-redux';
import { bindActionCreators } from 'redux';
import Snackbar from '../UI/Snackbar';
import type { Notifications } from '../../types/notifications';
import type { StateType } from '../../reducers/types';
import * as NotificationsActions from '../../actions/notifications';

type NotificationProps = {
  notifications: Notifications,
  closeNotification: string => void,
  destroyNotification: string => void
};

function NotificationsList({
  notifications,
  closeNotification,
  destroyNotification
}: NotificationProps) {
  const makeCloseNotification = k => () => {
    closeNotification(k);
    destroyNotification(k);
  };

  return (
    <>
      {Object.keys(notifications).map(k => (
        <Snackbar
          key={notifications[k].id}
          setClosed={makeCloseNotification(k)}
          message={notifications[k].message}
          variant={notifications[k].variant}
          isOpen={notifications[k].shown}
          duration={notifications[k].duration || 3000}
        />
      ))}
    </>
  );
}

// $FlowFixMe: flow disabled for this line
export default connect(
  (state: StateType) => ({
    notifications: state.notifications
  }),
  dispatch =>
    bindActionCreators(
      {
        closeNotification: NotificationsActions.closeNotification,
        destroyNotification: NotificationsActions.destroyNotification
      },
      dispatch
    )
)(NotificationsList);
