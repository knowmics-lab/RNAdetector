// @flow
import uniqid from 'uniqid';
import { is } from 'electron-util';
import { ipcRenderer } from 'electron';

const registeredCallbacks = new Map();

if (is.renderer) {
  ipcRenderer.on(
    'upload-message',
    (
      event,
      {
        id,
        isDone,
        isProgress,
        isError,
        error,
        percentage,
        bytesUploaded,
        bytesTotal
      }
    ) => {
      if (registeredCallbacks.has(id)) {
        // $FlowFixMe: undefined is checked with the previous line
        const { resolve, reject, onProgress } = registeredCallbacks.get(id);
        if (isDone) {
          registeredCallbacks.delete(id);
          resolve();
        } else if (isError) {
          reject(new Error(error));
        } else if (isProgress) {
          onProgress(percentage, bytesUploaded, bytesTotal);
        }
      }
    }
  );
}

export default {
  async upload(
    endpoint: string,
    filePath: string,
    fileName: string,
    fileType: string,
    onProgress: ?(number, number, number) => void
  ): Promise<void> {
    return new Promise((resolve, reject) => {
      const id = uniqid();
      registeredCallbacks.set(id, {
        resolve,
        reject,
        onProgress: (a, b, c) => {
          if (onProgress) onProgress(a, b, c);
        }
      });
      ipcRenderer.send('upload-file', {
        id,
        filePath,
        fileName,
        fileType,
        endpoint
      });
    });
  }
};
