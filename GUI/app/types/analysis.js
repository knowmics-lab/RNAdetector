/* eslint-disable flowtype/space-after-type-colon */
// @flow

import type { Job } from './jobs';

export type AnalysisFileTypes = 'fastq' | 'bam' | 'sam';

export type TrimGaloreConfig = {|
  enable?: boolean,
  quality?: number,
  length?: number,
  hardTrim?: boolean
|};

export type LongRNAAnalysisConfig = {|
  paired?: boolean,
  firstInputFile: string,
  secondInputFile: ?string,
  inputType: AnalysisFileTypes,
  convertBam?: boolean,
  trimGalore?: TrimGaloreConfig,
  algorithm?: 'salmon' | 'tophat' | 'hisat2',
  countingAlgorithm?: 'htseq' | 'feature-counts' | 'salmon',
  genome?: string,
  transcriptome?: string,
  annotation?: string,
  threads?: number
|};

export type SmallRNAAnalysisConfig = {|
  paired?: boolean,
  firstInputFile: string,
  secondInputFile: string,
  inputType: 'fastq' | 'bam' | 'sam',
  convertBam?: boolean,
  trimGalore?: TrimGaloreConfig,
  algorithm?: 'salmon' | 'tophat' | 'hisat2',
  countingAlgorithm?: 'htseq' | 'feature-counts' | 'salmon',
  genome?: string,
  transcriptome?: string,
  annotation?: string,
  threads?: number
|};

export type CircRNAAnalysisConfig = {|
  paired?: boolean,
  firstInputFile: string,
  secondInputFile: string,
  inputType: 'fastq' | 'bam' | 'sam',
  convertBam?: boolean,
  trimGalore?: TrimGaloreConfig,
  ciriQuant?: boolean,
  genome?: string,
  annotation?: string,
  bedAnnotation?: string,
  threads?: number,
  useFastqPair?: boolean,
  ciriSpanningDistance?: number,
  ciriUseVersion1?: boolean
|};

export type ContrastType = {|
  case: string,
  control: string
|};

export type SampleTypes = 'gene' | 'transcript';

export type DiffExpParameters = {|
  pcut?: number,
  log_offset?: number,
  when_apply_filter?: 'prenorm' | 'postnorm',
  norm?: 'edger' | 'deseq',
  norm_args?: {|
    method?: 'TMM' | 'TMMwsp' | 'RLE' | 'upperquartile' | 'none',
    locfunc?: 'median' | 'shorth'
  |},
  stats?: ('limma' | 'deseq' | 'edger')[],
  stats_args?: {|
    deseq?: {|
      fitType?: 'parametric' | 'local' | 'mean'
    |},
    edger?: {|
      main_method?: 'classic' | 'glm',
      rowsum_filter?: number,
      trend?: 'movingave' | 'loess' | 'none',
      tag_method?: 'grid' | 'optimize',
      glm_method?: 'CoxReid' | 'Pearson' | 'deviance',
      trend_method?: 'auto' | 'bin.spline' | 'bin.loess' | 'power' | 'spline'
    |},
    limma?: {|
      normalize_method?: 'none' | 'scale' | 'quantile' | 'cyclicloess'
    |}
  |},
  filters?: {|
    length?: {|
      length?: ?number
    |},
    avg_reads?: ?{|
      average_per_bp: number,
      quantile: number
    |},
    expression?: ?{|
      median?: boolean,
      mean?: boolean,
      quantile?: ?number,
      known?: ?(string[])
    |},
    presence?: ?{|
      frac: number,
      min_count: number,
      per_condition?: boolean
    |}
  |},
  adjust_method:
    | 'qvalue'
    | 'holm'
    | 'hochberg'
    | 'hommel'
    | 'bonferroni'
    | 'BH'
    | 'BY'
    | 'none',
  meta_p_method:
    | 'simes'
    | 'bonferroni'
    | 'minp'
    | 'maxp'
    | 'weight'
    | 'pandora'
    | 'dperm.min'
    | 'dperm.max'
    | 'dperm.weight'
    | 'fisher'
    | 'fperm'
    | 'whitlock'
    | 'none',
  fig_formats: ('png' | 'jpg' | 'tiff' | 'bmp' | 'pdf')[],
  num_cores: number
|};

export type DiffExpAnalysisConfig = {|
  source_sample_group: Job | number,
  sample_type: SampleTypes,
  condition_variables: string[],
  contrasts: ContrastType[],
  parameters: DiffExpParameters
|};
