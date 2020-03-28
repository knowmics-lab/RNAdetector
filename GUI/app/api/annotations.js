/* eslint-disable camelcase */
// @flow
import Connector from './connector';
import type {
  Annotation,
  AnnotationsCollection,
  AnnotationType
} from '../types/annotations';
import type { ResponseType, SimpleMapType, SortingSpec } from '../types/common';
import type { Job } from '../types/jobs';

export default {
  async create(
    name: string,
    type: AnnotationType,
    file: string,
    map_file: ?string = null
  ): Promise<ResponseType<Job>> {
    const result = await Connector.callPost('jobs', {
      name: `Create annotation ${name}`,
      type: 'annotation_upload_job_type',
      parameters: {
        name,
        type,
        file,
        map_file
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
  },
  async fetchAllByType(type: AnnotationType): Promise<SimpleMapType<string>> {
    const result = await Connector.callGet(`annotations`, {
      per_page: 0,
      type
    });
    const { data } = result.data;
    return Object.fromEntries(
      // $FlowFixMe: data is of type Annotation[]
      data.map(({ name }) => [name, name])
    );
  }
};
