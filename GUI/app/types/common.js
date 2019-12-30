// @flow

export type AxiosHeaders = {|
  headers: { [string]: string }
|};

export type ModifiableStateType = {
  +saving: boolean
};

export type MetaResponseType = {
  current_page: number,
  last_page: number,
  per_page: number,
  from: number,
  to: number,
  total: number,
  sorting?: SortingSpec
};

export type StatePaginationType = {|
  +current_page: ?number,
  +last_page: ?number,
  +per_page: number,
  +total: ?number,
  +sorting?: SortingSpec,
  +fetching: boolean
|};

export type LoadedCollectionMeta = {|
  +fetching: boolean
|};

export type MapType = { [string]: string | number | boolean | MapType };

export type SortingDirection = 'asc' | 'desc';

export type SortingSpec = { [string]: SortingDirection };
