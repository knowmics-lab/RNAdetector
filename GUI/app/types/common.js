// @flow

export type FileFilter = {
  name: string,
  extensions: string[]
};

export type DialogOptions = {
  title?: string,
  buttonLabel?: string,
  filters?: FileFilter[],
  message?: string,
  properties?: Array<'openFile' | 'openDirectory' | 'multiSelections'>
};

export type AxiosHeaders = {|
  headers: { [string]: mixed }
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

export type SimpleMapType<T> = { [string]: T };

export type RecursiveMapType<T> = { [string]: T | RecursiveMapType<T> };

export type MapType = {
  [string]: string | number | boolean | MapType | MapType[]
};

export type SortingDirection = 'asc' | 'desc';

export type SortingSpec = { [string]: SortingDirection };

export type ResponseType<T> = {
  validationErrors?: RecursiveMapType<string>,
  data?: T
};

export type UploadProgressFunction = (number, number, number) => void;
