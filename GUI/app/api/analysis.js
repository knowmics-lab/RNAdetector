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
  name,
  type,
  parameters
): Promise<ResponseType<Job>> {
  const result = await Connector.callPost('jobs', {
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
  async createLongRNA(
    name: string,
    parameters: LongRNAAnalysisConfig
  ): Promise<ResponseType<Job>> {
    return realJobSubmit(name, 'long_rna_job_type', parameters);
  },
  async createSmallRNA(
    name: string,
    parameters: SmallRNAAnalysisConfig
  ): Promise<ResponseType<Job>> {
    return realJobSubmit(name, 'small_rna_job_type', parameters);
  },
  async createCircRNA(
    name: string,
    parameters: CircRNAAnalysisConfig
  ): Promise<ResponseType<Job>> {
    return realJobSubmit(name, 'circ_rna_job_type', parameters);
  }
};
