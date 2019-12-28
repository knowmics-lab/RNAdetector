// @flow
import axios from 'axios';
import Settings from './settings';
import type { Job, JobsCollection, JobTypesCollection } from '../types/jobs';

export default {
  async processDeletedList(deleted: number[]): Promise<number[]> {
    if (deleted.length === 0) return deleted;
    const deletedPromises = deleted.map(
      id =>
        new Promise(resolve => {
          this.fetchJobById(id)
            .then(() => resolve(true))
            .catch(() => resolve(false));
        })
    );
    const res = await Promise.all(deletedPromises);
    return deleted.filter((ignore, idx) => res[idx]);
  },
  async deleteJob(jobId: number): Promise<void> {
    await axios.delete(`${Settings.getApiUrl()}jobs/${jobId}`, {
      ...Settings.getAxiosHeaders()
    });
  },
  async submitJob(jobId: number): Promise<Job> {
    const result = await axios.get(
      `${Settings.getApiUrl()}jobs/${jobId}/submit`,
      {
        ...Settings.getAxiosHeaders()
      }
    );
    const { data, links } = result.data;
    return {
      ...data,
      links
    };
  },
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
      data: data.map(x => ({
        ...x,
        links: {
          self: x.self,
          upload: x.upload,
          submit: x.submit
        }
      })),
      meta
    };
  },
  async fetchTypes(): Promise<JobTypesCollection> {
    const result = await axios.get(`${Settings.getApiUrl()}job-types`, {
      ...Settings.getAxiosHeaders()
    });
    const { data } = result.data;
    return data;
  }
};
