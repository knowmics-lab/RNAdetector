// @flow
import type {
  MetaResponseType,
  SimpleMapType,
  StatePaginationType
} from './common';

export type IndexingAlgorithm = 'bwa' | 'hisat' | 'salmon' | 'star';

export type Reference = {
  id: number,
  name: string,
  available_for: SimpleMapType<boolean>,
  created_at: string,
  created_at_diff: string,
  updated_at: string,
  updated_at_diff: string,
  links?: {
    self: string
  }
};

export type ReferencesCollection = {
  data: Reference[],
  meta: MetaResponseType
};

export type ReferencesListType = {|
  +refreshAll: boolean,
  +refreshPages: number[],
  +state: StatePaginationType,
  +pages: { +[number]: Reference[] }
|};

export type LoadedReferences = {|
  fetching: boolean,
  +items: { +[number]: Reference }
|};

export type ReferencesStateType = {|
  +referencesList: ReferencesListType,
  +references: LoadedReferences
|};
