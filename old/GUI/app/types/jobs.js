// @flow
import type { MapType, MetaResponseType, StatePaginationType } from './common';

export type JobPathType = {
  path: string,
  url: string
};

export type MetadataType = {
  field: string,
  type: string,
  content: string[]
};

export type ConfirmationOutputType = {
  type: 'confirmation',
  done: boolean
};

export type AnalysisOutputType = {
  type: 'analysis',
  outputFile: JobPathType
};

export type HarmonizedAnalysisOutputType = {
  type: 'analysis-harmonized',
  outputFile: JobPathType,
  harmonizedFile: JobPathType
};

export type HarmonizedTranscriptsAnalysisOutputType = {
  type: 'analysis-harmonized-transcripts',
  outputFile: JobPathType,
  harmonizedFile: JobPathType,
  harmonizedTranscriptsFile: JobPathType
};

export type HarmonizedDescriptionAnalysisOutputType = {
  type: 'analysis-harmonized-description',
  outputFile: JobPathType,
  harmonizedFile: JobPathType,
  description: JobPathType,
  metadata?: MetadataType[]
};

export type HarmonizedTranscriptsDescriptionAnalysisOutputType = {
  type: 'analysis-harmonized-transcripts-description',
  outputFile: JobPathType,
  harmonizedFile: JobPathType,
  harmonizedTranscriptsFile: JobPathType,
  description: JobPathType,
  metadata?: MetadataType[]
};

export type AnalysisReportOutputType = {
  type: 'analysis-report',
  outputFile: JobPathType,
  reportFile: JobPathType
};

export type JobOutputType =
  | ConfirmationOutputType
  | AnalysisOutputType
  | AnalysisReportOutputType
  | HarmonizedAnalysisOutputType
  | HarmonizedTranscriptsAnalysisOutputType
  | HarmonizedDescriptionAnalysisOutputType
  | HarmonizedTranscriptsDescriptionAnalysisOutputType;

export type Job = {
  id: number,
  sample_code?: string,
  name: string,
  type: string,
  readable_type: string,
  status: 'ready' | 'queued' | 'processing' | 'completed' | 'failed',
  parameters?: MapType,
  output?: JobOutputType,
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
  fetching: boolean,
  submitting: number[],
  deleting: number[],
  +items: { +[number]: Job }
|};

export type JobsStateType = {|
  +jobsList: JobsListType,
  +jobs: LoadedJobs
|};
