import type { Dispatch as ReduxDispatch, Store as ReduxStore } from 'redux';
import type { ConfigObjectType, JobsCollectionItem } from '../api';

export type counterStateType = {
  +counter: number
};

export type settingsStateType = {
  +settings: {
    +state: {
      +saving: boolean,
      +saved: boolean,
      +error: boolean,
      +message: string
    }
  } & ConfigObjectType
};

export type jobsListType = {
  +jobsList: {|
    +state: {
      +current_page: number,
      +last_page: number,
      +per_page: number,
      +total: number,
      +isFetching: boolean,
      +isError: boolean,
      +errorMessage: string
    },
    +pages: { +[number]: JobsCollectionItem[] }
  |}
};

export type Action = {
  +type: string,
  payload?: {}
};

export type StateType = counterStateType & settingsStateType & jobsListType;

export type GetState = () => StateType;

export type Dispatch = ReduxDispatch<Action>;

export type Store = ReduxStore<GetState, Action>;
