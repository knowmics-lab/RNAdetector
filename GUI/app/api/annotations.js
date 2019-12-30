/* eslint-disable camelcase */
// @flow
import axios from 'axios';
import Settings from './settings';
import type { Annotation, AnnotationsCollection } from '../types/annotations';
import type { SortingSpec } from '../types/common';

export default {
  async delete(id: number): Promise<void> {
    await axios.delete(`${Settings.getApiUrl()}annotations/${id}`, {
      ...Settings.getAxiosHeaders()
    });
  },
  async fetchById(id: number): Promise<Annotation> {
    const result = await axios.get(`${Settings.getApiUrl()}annotations/${id}`, {
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
  ): Promise<AnnotationsCollection> {
    const order = Object.keys(sorting);
    const order_direction = Object.values(sorting);
    const result = await axios.get(`${Settings.getApiUrl()}annotations`, {
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
