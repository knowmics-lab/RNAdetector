// @flow
import axios from 'axios';
import Settings from './settings';
import type { Job, JobsCollection } from '../types/jobs';

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
