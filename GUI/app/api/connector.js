/* eslint-disable no-restricted-syntax,no-plusplus */
import axios from 'axios';
import { set } from 'lodash';
import Settings from './settings';
import type { MapType, ResponseType } from '../types/common';

const dotParser = <T>(obj: T, suppress = []) => {
  const rxp = new RegExp(suppress.join('|'));
  const res = {};
  for (const k of Object.keys(obj)) {
    const sanitized = k.replace(rxp, '').replace(/^\./, '');
    const val = Array.isArray(obj[k]) && obj[k].length ? obj[k][0] : obj[k];
    if (sanitized) {
      set(res, sanitized, val);
    }
  }
  return res;
};

const parseErrorResponse = <T>(response): ResponseType<T> => {
  if (response.status && response.status === 422) {
    const validationErrors = dotParser(response.data.errors, ['^parameters']);
    return {
      validationErrors
    };
  }
  if (response.data && response.data.message) {
    throw new Error(response.data.message);
  }
  throw new Error('Unknown error');
};

export default {
  getEndpointUrl(endpoint: string): string {
    return `${Settings.getApiUrl()}${endpoint.replace(/^\//gm, '')}`;
  },
  async callGet<T>(
    endpoint: string,
    params: MapType = {}
  ): Promise<ResponseType<T>> {
    const url = this.getEndpointUrl(endpoint);
    try {
      const { data } = await axios.get(url, {
        params,
        ...Settings.getAxiosHeaders()
      });
      return {
        data
      };
    } catch (e) {
      if (e.response) {
        const { response } = e;
        return parseErrorResponse(response);
      }
      throw e;
    }
  },
  async callPost<T>(
    endpoint: string,
    params: MapType = {}
  ): Promise<ResponseType<T>> {
    const url = this.getEndpointUrl(endpoint);
    try {
      const { data } = await axios.post(url, params, {
        ...Settings.getAxiosHeaders()
      });
      return {
        data
      };
    } catch (e) {
      if (e.response) {
        const { response } = e;
        return parseErrorResponse(response);
      }
      throw e;
    }
  },
  async callPatch<T>(
    endpoint: string,
    params: MapType = {}
  ): Promise<ResponseType<T>> {
    const url = this.getEndpointUrl(endpoint);
    try {
      const { data } = await axios.patch(url, params, {
        ...Settings.getAxiosHeaders()
      });
      return {
        data
      };
    } catch (e) {
      if (e.response) {
        const { response } = e;
        return parseErrorResponse(response);
      }
      throw e;
    }
  },
  async callDelete<T>(
    endpoint: string,
    params: MapType = {}
  ): Promise<ResponseType<T>> {
    const url = this.getEndpointUrl(endpoint);
    try {
      const { data } = await axios.delete(url, {
        params,
        ...Settings.getAxiosHeaders()
      });
      return {
        data
      };
    } catch (e) {
      if (e.response) {
        const { response } = e;
        return parseErrorResponse(response);
      }
      throw e;
    }
  }
};
