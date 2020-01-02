/* eslint-disable camelcase */
// @flow
import path from 'path';
import fs from 'fs-extra';
import { api as electron } from 'electron-util';
import Settings from './settings';
import type { Job, JobsCollection, JobTypesCollection } from '../types/jobs';
import type { SortingSpec } from '../types/common';
import Connector from './connector';
import Downloader from './downloader';

export default {
  getUploadUrl(job: number | Job): string {
    let jobId;
    if (typeof job === 'object') {
      jobId = job.id;
    } else {
      jobId = job;
    }
    return Connector.getEndpointUrl(`jobs/${jobId}/upload`);
  },
  async download(
    jobId: number,
    onStart?: () => void,
    onCompleted?: () => void
  ): Promise<void> {
    const job = await this.fetchJobById(jobId);
    if (job.output && job.output.outputFile && job.output.outputFile.path) {
      const { path: outputPath } = job.output.outputFile;
      if (typeof outputPath === 'string') {
        const outputUrl = Settings.getPublicUrl(outputPath);
        const outputFilename = path.basename(outputPath);
        Downloader.downloadUrl(outputUrl, outputFilename, onStart, onCompleted);
      }
    } else {
      throw new Error('Unable to find output path');
    }
  },
  async openLocalFolder(jobId: number): Promise<void> {
    if (!Settings.isLocal())
      throw new Error('Unable to open folder when in remote mode');
    const job = await this.fetchJobById(jobId);
    if (job.output && job.output.outputFile && job.output.outputFile.path) {
      const { path: outputPath } = job.output.outputFile;
      if (typeof outputPath === 'string') {
        const jobFolder = path.dirname(
          Settings.getLocalPath(`/public/${outputPath}`)
        );
        if (!(await fs.pathExists(jobFolder)))
          throw new Error('Unable to find output path');
        if (!electron.shell.openItem(jobFolder)) {
          throw new Error('Unable to open output folder');
        }
      }
    } else {
      throw new Error('Unable to find output path');
    }
  },
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
    await Connector.callDelete(`jobs/${jobId}`);
  },
  async submitJob(jobId: number): Promise<Job> {
    const result = await Connector.callGet(`jobs/${jobId}/submit`);
    const { data, links } = result.data;
    return {
      ...data,
      links
    };
  },
  async fetchJobById(jobId: number): Promise<Job> {
    const result = await Connector.callGet(`jobs/${jobId}`);
    const { data, links } = result.data;
    return {
      ...data,
      links
    };
  },
  async fetchJobs(
    per_page: number = 15,
    sorting: SortingSpec = { created_at: 'desc' },
    page: number = 1
  ): Promise<JobsCollection> {
    const order = Object.keys(sorting);
    const order_direction = Object.values(sorting);
    const result = await Connector.callGet(`jobs`, {
      page,
      per_page,
      order,
      order_direction
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
      meta: {
        ...meta,
        sorting
      }
    };
  },
  async fetchTypes(): Promise<JobTypesCollection> {
    const result = await Connector.callGet(`job-types`);
    const { data } = result.data;
    return data;
  }
};
