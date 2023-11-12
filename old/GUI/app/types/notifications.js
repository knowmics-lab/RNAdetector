// @flow

export type TypeOfNotification = 'success' | 'warning' | 'error' | 'info';

export type Notification = {
  message: string,
  variant: TypeOfNotification,
  duration: ?number
};

export type PushedNotification = Notification & { id: string, shown: boolean };

export type Notifications = { [string]: PushedNotification };

export type NotificationsState = {
  +notifications: Notifications
};

export type PushNotificationFunction = (string, ?TypeOfNotification) => void;
