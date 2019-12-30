// @flow
import type { MetaResponseType, StatePaginationType } from './common';

export type Annotation = {
  id: number,
  name: string,
  created_at: string,
  created_at_diff: string,
  updated_at: string,
  updated_at_diff: string,
  links?: {
    self: string
  }
};

export type AnnotationsCollection = {
  data: Annotation[],
  meta: MetaResponseType
};

export type AnnotationsListType = {|
  +refreshAll: boolean,
  +refreshPages: number[],
  +state: StatePaginationType,
  +pages: { +[number]: Annotation[] }
|};

export type LoadedAnnotations = {|
  fetching: boolean,
  +items: { +[number]: Annotation }
|};

export type AnnotationsStateType = {|
  +annotationsList: AnnotationsListType,
  +annotations: LoadedAnnotations
|};
