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
    code: string,
    name: string,
    jobs: (Job | number)[],
    descriptionFile?: ?string,
    de_novo: boolean = false
  ): Promise<ResponseType<Job>> {
    return realJobSubmit(code, name, 'samples_group_job_type', {
      jobs: jobs.map(j => (typeof j === 'object' ? j.id : j)),
      description: descriptionFile || undefined,
      de_novo
    });
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
