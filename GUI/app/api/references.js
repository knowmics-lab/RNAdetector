/* eslint-disable camelcase */
// @flow
import Connector from './connector';
import type {
  IndexingAlgorithm,
  Reference,
  ReferencesCollection
} from '../types/references';
import type { SortingSpec, ResponseType, SimpleMapType } from '../types/common';
import type { Job } from '../types/jobs';
import type { AnnotationType } from '../types/annotations';

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
    await Connector.callDelete(`references/${id}`);
  },
  async fetchById(id: number): Promise<Reference> {
    const result = await Connector.callGet(`references/${id}`);
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
    const result = await Connector.callGet(`references`, {
      page,
      per_page,
      order,
      order_direction
    });
    const { data, meta } = result.data;
    return {
      data,
      meta: {
        ...meta,
        sorting
      }
    };
  },
  async fetchAllByAlgorithm(
    algorithm: IndexingAlgorithm
  ): Promise<SimpleMapType<string>> {
    const result = await Connector.callGet(`references`, {
      per_page: 0,
      indexed_for: algorithm
    });
    const { data } = result.data;
    return Object.fromEntries(
      // $FlowFixMe: data is of type Reference[]
      data.map(({ name }) => [name, name])
    );
  }
};
