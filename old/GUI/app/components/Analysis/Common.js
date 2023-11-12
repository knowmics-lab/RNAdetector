import * as Api from '../../api';
import type { File } from '../UI/FileSelector2';
import type { Job } from '../../types/jobs';
import type { ResponseType } from '../../types/common';
import { RecursiveMapType } from '../../types/common';

export default {
  filterParamsByAlgorithm(params: *) {
    if (params.algorithm === 'salmon') {
      return Api.Utils.filterByKey(
        params,
        k => k !== 'annotation' && k !== 'genome'
      );
    }

    if (
      params.annotation !== 'salmon' &&
      params.countingAlgorithm !== 'salmon'
    ) {
      return Api.Utils.filterByKey(params, k => k !== 'transcriptome');
    }

    return params;
  },
  checkLength(
    files: [?File, ?File][],
    paired: boolean,
    pushNotification: (
      string,
      ?('success' | 'warning' | 'error' | 'info')
    ) => void
  ): number {
    const validLength = files.filter(
      f => !!f[0] && (!paired || (paired && !!f[1]))
    ).length;
    const firstLength = files.map(f => f[0]).filter(f => !!f).length;
    const secondLength = files.map(f => f[1]).filter(f => !!f).length;
    if (validLength < 1) {
      pushNotification('You should select at least one input file.', 'error');
      return -1;
    }
    if (
      paired &&
      (firstLength !== validLength || secondLength !== validLength)
    ) {
      pushNotification(
        'You must select the same number of mate input files.',
        'error'
      );
      return -1;
    }

    return validLength;
  },
  async uploadFile(job: Job, file: ?File, setState: (*) => void) {
    if (!file) return;
    Api.Upload.ui.uploadStart(setState, file.name);
    await Api.Upload.upload(
      job,
      file.path,
      file.name,
      file.type,
      Api.Upload.ui.makeOnProgress(setState)
    );
    Api.Upload.ui.uploadEnd(setState);
  },
  async createGroup(
    single: boolean,
    code: string,
    name: string,
    jobs: Job[],
    descriptionFile: ?File,
    pushNotification: (
      string,
      ?('success' | 'warning' | 'error' | 'info')
    ) => void,
    handleValidationError: (boolean, RecursiveMapType<string>) => void,
    setState: (*) => void
  ): Promise<?Job> {
    if (single) return null;
    const data: ResponseType<Job> = await Api.Analysis.createSampleGroup(
      code,
      name,
      jobs,
      descriptionFile ? descriptionFile.name : undefined
    );
    if (data.validationErrors) {
      pushNotification(
        'Errors occurred during validation of input parameters. Please review the form!',
        'warning'
      );
      handleValidationError(false, data.validationErrors);
      return undefined;
    }
    const { data: job } = data;
    if (descriptionFile) {
      await this.uploadFile(job, descriptionFile, setState);
    }
    pushNotification('Sample group created!');
    return job;
  }
};
