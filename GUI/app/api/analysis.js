/* eslint-disable camelcase */
// @flow

import type { ResponseType } from '../types/common';
import type { Job } from '../types/jobs';
import Connector from './connector';
import type {
  CircRNAAnalysisConfig,
  LongRNAAnalysisConfig,
  SmallRNAAnalysisConfig
} from '../types/analysis';

async function realJobSubmit(
  sample_code,
  name,
  type,
  parameters
): Promise<ResponseType<Job>> {
  const result = await Connector.callPost('jobs', {
    sample_code,
    name,
    type,
    parameters
  });
  if (result.validationErrors) {
    return result;
  }
  const { data, links } = result.data;
  return {
    data: {
      ...data,
      links
    }
  };
}

export default {
  async createSampleGroup(
    name: string,
    
  ): Promise<ResponseType<Job>> {

  },
  async createLongRNA(
    code: string,
    name: string,
    parameters: LongRNAAnalysisConfig
  ): Promise<ResponseType<Job>> {
    return realJobSubmit(code, name, 'long_rna_job_type', parameters);
  },
  async createSmallRNA(
    code: string,
    name: string,
    parameters: SmallRNAAnalysisConfig
  ): Promise<ResponseType<Job>> {
    return realJobSubmit(code, name, 'small_rna_job_type', parameters);
  },
  async createCircRNA(
    code: string,
    name: string,
    parameters: CircRNAAnalysisConfig
  ): Promise<ResponseType<Job>> {
    return realJobSubmit(code, name, 'circ_rna_job_type', parameters);
  }
};
