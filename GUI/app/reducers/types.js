import type { Dispatch as ReduxDispatch, Store as ReduxStore } from 'redux';
import type { ConfigObjectType } from '../api';

export type counterStateType = {
  +counter: number
};

export type settingsStateType = {
  +settings: {
    +state: {
      saving: boolean,
      saved: boolean,
      error: boolean,
      message: string
    }
  } & ConfigObjectType
};

export type Action = {
  +type: string,
  payload?: {}
};

export type GetState = () => counterStateType & settingsStateType;

export type Dispatch = ReduxDispatch<Action>;

export type Store = ReduxStore<GetState, Action>;
