// @flow
import Settings from './settings';
import axios from 'axios';

export type Job = {
  id: number,
  type: string,
  status: 'ready' | 'queued' | 'processing' | 'completed' | 'failed',
  parameters: { [string]: string | number | {} },
  output: { [string]: string | number | {} },
  log: string,
  created_at: string,
  updated_at: string,
  owner: *,
  links: {
    self: string,
    owner: string,
    upload: string,
    submit: string
  }
};

export type JobsCollectionItem = {
  id: string,
  type: string,
  status: 'ready' | 'queued' | 'processing' | 'completed' | 'failed',
  created_at: string,
  updated_at: string,
  owner: *,
  self: string,
  upload: string,
  submit: string
};

export type JobsCollection = {
  data: { [number]: JobsCollectionItem },
  meta: {
    current_page: number,
    last_page: number,
    per_page: number,
    from: number,
    to: number,
    total: number
  }
};

export default {
  async fetchJobById(jobId: number): Promise<Job> {
    const result = await axios.get(`${Settings.getApiUrl()}jobs/${jobId}`, {
      ...Settings.getAxiosHeaders()
    });
    const { data, links } = result.data;
    return {
      ...data,
      links
    };
  },
  async fetchJobs(
    per_page: number = 15,
    page: number = 1
  ): Promise<JobsCollection> {
    const result = await axios.get(`${Settings.getApiUrl()}jobs`, {
      params: {
        page,
        per_page
      },
      ...Settings.getAxiosHeaders()
    });
    const { data, meta } = result.data;
    return {
      data,
      meta
    };
  }
};
