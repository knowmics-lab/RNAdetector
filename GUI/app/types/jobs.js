// @flow
import type {
  LoadedCollectionMeta,
  MapType,
  MetaResponseType,
  StatePaginationType
} from './common';

export type Job = {
  id: number,
  name: string,
  type: string,
  readable_type: string,
  status: 'ready' | 'queued' | 'processing' | 'completed' | 'failed',
  parameters?: MapType,
  output?: MapType,
  log?: string,
  created_at: string,
  created_at_diff: string,
  updated_at: string,
  updated_at_diff: string,
  owner: *,
  links: {
    self: string,
    owner?: string,
    upload: string,
    submit: string
  }
};

export type JobsCollectionItem = Job;

export type JobsCollection = {
  data: JobsCollectionItem[],
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

export type LoadedJobs = {|
  +meta: LoadedCollectionMeta,
  +items: { +[number]: Job }
|};

export type JobsStateType = {|
  +jobsList: JobsListType,
  +jobs: LoadedJobs
|};
