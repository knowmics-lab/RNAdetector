// @flow

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
  secondInputFile: string,
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
  genome?: string,
  annotation?: string,
  threads?: number,
  ciriSpanningDistance?: number
|};
