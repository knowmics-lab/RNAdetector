/* eslint-disable camelcase */
// @flow
import Connector from './connector';
import type {
  IndexingAlgorithm,
  Reference,
  ReferencesCollection
} from '../types/references';
import type {
  SortingSpec,
  ResponseType,
  SimpleMapType,
  MapType,
  RecursiveMapType
} from '../types/common';
import type { Job } from '../types/jobs';
import type { Package } from '../types/local';

export default {
  async listPackages(): Promise<Package[]> {
    const result = await Connector.callGet(`references/packages`);
    return result.data.packages;
  },
  async install(names: string[]): Promise<ResponseType<Job>> {
    const result = await Connector.callPost('jobs', {
      name:
        names.length === 1 ? `Install package ${names[0]}` : 'Install packages',
      type: 'install_packages_job_type',
      parameters: {
        names
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
  async create(
    name: string,
    fastaFile: string,
    availableFor: string[],
    map_file: ?string = null,
    customArguments: ?SimpleMapType<string> = undefined
  ): Promise<ResponseType<Job>> {
    const parameters: RecursiveMapType<*> = {
      name,
      fastaFile,
      index: Object.fromEntries(
        ['bwa', 'hisat', 'salmon', 'star'].map(v => [
          v,
          availableFor.includes(v)
        ])
      ),
      map_file
    };
    if (customArguments && customArguments !== {}) {
      parameters.custom_arguments = customArguments;
    }
    const result = await Connector.callPost('jobs', {
      name: `Create and index reference ${name}`,
      type: 'reference_upload_job_type',
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
