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

//export type PaginatedJobsCollection

export type JobsCollection = { [number]: JobsCollectionItem };

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
  async fetchJobs(): Promise<JobsCollection> {

  }
};
