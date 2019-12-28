// @flow
import type { Action } from './types';
import * as NotificationsActions from '../actions/notifications';
import type { Notifications } from '../types/notifications';

export default function notifications(
  state: Notifications,
  action: Action
): Notifications {
  const oldState = state || {};
  switch (action.type) {
    case NotificationsActions.NOTIFICATION_PUSH:
      return {
        ...state,
        [action.payload.id]: action.payload
      };
    case NotificationsActions.NOTIFICATION_CLOSE:
      return {
        ...state,
        [action.payload]: {
          ...state[action.payload],
          shown: false
        }
      };
    case NotificationsActions.NOTIFICATION_DESTROY:
      return Object.keys(state).reduce((object, key) => {
        if (key !== action.payload) {
          // eslint-disable-next-line no-param-reassign
          object[key] = state[key];
        }
        return object;
      }, {});
    default:
      return oldState;
  }
}
