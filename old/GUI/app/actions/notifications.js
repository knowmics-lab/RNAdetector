// @flow
import uniqid from 'uniqid';
import type { Notification } from '../types/notifications';
import type { Action } from '../reducers/types';

export const NOTIFICATION_PUSH = 'NOTIFICATIONS--PUSH';
export const NOTIFICATION_CLOSE = 'NOTIFICATIONS--CLOSE';
export const NOTIFICATION_DESTROY = 'NOTIFICATIONS--DESTROY';

export function pushNotificationSimple(
  message: string,
  variant: 'success' | 'warning' | 'error' | 'info' = 'success'
): Action {
  return pushNotification({
    message,
    variant,
    duration: null
  });
}

export function pushNotification(payload: Notification): Action {
  return {
    type: NOTIFICATION_PUSH,
    payload: {
      ...payload,
      id: uniqid(),
      shown: true
    }
  };
}

export function closeNotification(payload: string): Action {
  return {
    type: NOTIFICATION_CLOSE,
    payload
  };
}

export function destroyNotification(payload: string): Action {
  return {
    type: NOTIFICATION_DESTROY,
    payload
  };
}
