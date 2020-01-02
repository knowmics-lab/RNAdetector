/* eslint-disable camelcase */
// @flow
import axios from 'axios';
import Settings from './settings';
import Connector from './connector';
import type { Reference, ReferencesCollection } from '../types/references';
import type { SortingSpec, ResponseType } from '../types/common';
import type { Job } from '../types/jobs';

export default {
  async create(
    name: string,
    fastaFile: string,
    availableFor: string[]
  ): Promise<ResponseType<Job>> {
    const result = await Connector.callPost('jobs', {
      name: `Create and index reference ${name}`,
      type: 'reference_upload_job_type',
      parameters: {
        name,
        fastaFile,
        index: Object.fromEntries(
          ['bwa', 'tophat', 'hisat', 'salmon'].map(v => [
            v,
            availableFor.includes(v)
          ])
        )
      }
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
  },
  async delete(id: number): Promise<void> {
    await axios.delete(`${Settings.getApiUrl()}references/${id}`, {
      ...Settings.getAxiosHeaders()
    });
  },
  async fetchById(id: number): Promise<Reference> {
    const result = await axios.get(`${Settings.getApiUrl()}references/${id}`, {
      ...Settings.getAxiosHeaders()
    });
    const { data, links } = result.data;
    return {
      ...data,
      links
    };
  },
  async fetch(
    per_page: number = 15,
    sorting: SortingSpec = { created_at: 'desc' },
    page: number = 1
  ): Promise<ReferencesCollection> {
    const order = Object.keys(sorting);
    const order_direction = Object.values(sorting);
    const result = await axios.get(`${Settings.getApiUrl()}references`, {
      params: {
        page,
        per_page,
        order,
        order_direction
      },
      ...Settings.getAxiosHeaders()
    });
    const { data, meta } = result.data;
    return {
      data,
      meta: {
        ...meta,
        sorting
      }
    };
  }
};
