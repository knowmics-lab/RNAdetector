/* eslint-disable camelcase */
// @flow
import Connector from './connector';
import type { Annotation, AnnotationsCollection } from '../types/annotations';
import type { ResponseType, SortingSpec } from '../types/common';
import type { Job } from '../types/jobs';

export default {
  async create(
    name: string,
    type: 'gtf' | 'bed',
    file: string
  ): Promise<ResponseType<Job>> {
    const result = await Connector.callPost('jobs', {
      name: `Create annotation ${name}`,
      type: 'annotation_upload_job_type',
      parameters: {
        name,
        type,
        file
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
    await Connector.callDelete(`annotations/${id}`);
  },
  async fetchById(id: number): Promise<Annotation> {
    const result = await Connector.callGet(`annotations/${id}`);
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
  ): Promise<AnnotationsCollection> {
    const order = Object.keys(sorting);
    const order_direction = Object.values(sorting);
    const result = await Connector.callGet(`annotations`, {
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
  }
};
