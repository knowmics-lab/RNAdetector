import type { Dispatch as ReduxDispatch, Store as ReduxStore } from 'redux';
import type { JobsStateType } from '../types/jobs';
import type { SettingsStateType } from '../types/settings';
import type { NotificationsState } from '../types/notifications';

export type counterStateType = {
  +counter: number
};

export type Action = {
  +type: string,
  payload?: {}
};

export type StateType = counterStateType &
  SettingsStateType & { jobs: JobsStateType } & NotificationsState;

export type GetState = () => StateType;

export type Dispatch = ReduxDispatch<Action>;

export type Store = ReduxStore<GetState, Action>;
