// @flow

export type Notification = {
  message: string,
  variant: 'success' | 'warning' | 'error' | 'info',
  duration: ?number
};

export type PushedNotification = Notification & { id: string, shown: boolean };

export type Notifications = { [string]: PushedNotification };

export type NotificationsState = {
  +notifications: Notifications
};
