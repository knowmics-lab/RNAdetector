// @flow

import Tus from '@uppy/tus';
import uniqid from 'uniqid';
import Uppy from '@uppy/core';
import * as Api from './index';

type NotificationFunction = (string, ?string) => void;

const instances = new Map();
const settings = new Map();

export default {
  initUppyInstance(
    allowedFileTypes: ?(string[]) = null,
    minNumberOfFiles: number = 1,
    maxNumberOfFiles: number = 1
  ) {
    const id = uniqid();
    instances.set(
      id,
      Uppy({
        restrictions: {
          maxNumberOfFiles,
          allowedFileTypes
        },
        allowMultipleUploads: false,
        autoProceed: false
      })
    );
    settings.set(id, [minNumberOfFiles, maxNumberOfFiles]);
    return id;
  },
  getInstance(id: string) {
    return instances.get(id);
  },
  getSettings(id: string) {
    return settings.get(id);
  },
  configureEndpoint(id: string, url: string) {
    const uppy = this.getInstance(id);
    if (uppy.getPlugin('Tus') !== null) {
      uppy.removePlugin(uppy.getPlugin('Tus'));
    }
    uppy.use(Tus, {
      endpoint: url,
      headers: {
        ...Api.Settings.getAuthHeaders()
      }
    });
  },
  isValid(id: string, notify: NotificationFunction) {
    const uppy = this.getInstance(id);
    const files = uppy.getFiles();
    const [min, max] = this.getSettings(id);
    if (files.length < min) {
      notify(
        `You must select at least ${min} file${min > 1 ? 's' : ''}.`,
        'error'
      );
    }
    if (files.length > max) {
      notify(
        `You must select at most ${max} file${max > 1 ? 's' : ''}.`,
        'error'
      );
    }
    return min <= files.length && files.length <= max;
  },
  getFilename(id: string, index: number): ?string {
    const uppy = this.getInstance(id);
    const file = uppy.getFiles()[index];
    if (file) {
      return file.data.name;
    }
    return null;
  },
  checkUppyResult(
    id: string,
    uploadResult: *,
    index: number,
    notify: NotificationFunction
  ) {
    if (
      uploadResult.successful &&
      Array.isArray(uploadResult.successful) &&
      uploadResult.successful.length > index
    ) {
      if (uploadResult.successful[index].name === this.getFilename(id, index)) {
        return true;
      }
      notify(
        `Error during upload: uploaded file ${index +
          1} is different from the selected one.`,
        'error'
      );
    } else {
      notify(`Unknown error during upload of the ${index + 1} file!`, 'error');
    }
    return false;
  },
  async upload(
    id: string,
    url: string,
    notify: NotificationFunction
  ): Promise<boolean> {
    const uppy = this.getInstance(id);
    this.configureEndpoint(id, url);
    const uploadResult = await uppy.upload();
    const [min] = this.getSettings(id);
    let res = true;
    for (let i = 0; i < min; i += 1) {
      res = res && this.checkUppyResult(id, uploadResult, i, notify);
    }
    return res;
  },
  clearInstance(id: string) {
    const instance = this.getInstance(id);
    instance.close();
    instances.delete(id);
    settings.delete(id);
  }
};
