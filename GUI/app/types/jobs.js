// @flow
import type { MetaResponseType, StatePaginationType } from './common';

export type Job = {
  id: number,
  type: string,
  status: 'ready' | 'queued' | 'processing' | 'completed' | 'failed',
  parameters: { [string]: string | number | {} },
  output: { [string]: string | number | {} },
  log: string,
  created_at: string,
  updated_at: string,
  owner: *,
  links: {
    self: string,
    owner: string,
    upload: string,
    submit: string
  }
};

export type JobsCollectionItem = {
  id: string,
  type: string,
  status: 'ready' | 'queued' | 'processing' | 'completed' | 'failed',
  created_at: string,
  updated_at: string,
  owner: *,
  self: string,
  upload: string,
  submit: string
};

export type JobsCollection = {
  data: { [number]: JobsCollectionItem },
  meta: MetaResponseType
};

export type JobTypes = {
  id: string,
  description: string
};

export type JobTypesCollection = { [string]: JobTypes };

export type JobsListType = {|
  +refreshAll: boolean,
  +refreshPages: number[],
  +state: StatePaginationType,
  +pages: { +[number]: JobsCollectionItem[] }
|};

export type JobsStateType = {|
  +jobsList: JobsListType
|};
