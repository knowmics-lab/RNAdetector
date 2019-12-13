import type { Dispatch as ReduxDispatch, Store as ReduxStore } from 'redux';

export type counterStateType = {
  +counter: number
};

export type settingsStateType = {
  +settings: {
    +local: boolean,
    +webserviceUrl: string,
    +jobsPath: string
  }
};

export type Action = {
  +type: string,
  data?: {}
};

export type GetState = () => counterStateType & settingsStateType;

export type Dispatch = ReduxDispatch<Action>;

export type Store = ReduxStore<GetState, Action>;
