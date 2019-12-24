// @flow

export type AxiosHeaders = {|
  headers: { [string]: string }
|};

export type ModifiableStateType = {
  +saving: boolean,
  +saved: boolean,
  +error: boolean,
  +message: string
};

export type MetaResponseType = {
  current_page: number,
  last_page: number,
  per_page: number,
  from: number,
  to: number,
  total: number
};

export type StatePaginationType = {
  +current_page: ?number,
  +last_page: ?number,
  +per_page: number,
  +total: ?number,
  +fetching: boolean,
  +fetched: boolean,
  +error: ?string
};

export type LoadedCollectionMeta = {
  +fetching: boolean,
  +fetched: boolean,
  +error: ?string
};

export type MapType = { [string]: string | number | boolean | MapType };
