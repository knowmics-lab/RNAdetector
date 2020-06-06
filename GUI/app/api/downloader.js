// @flow
import { ipcRenderer } from 'electron';
import { is } from 'electron-util';
import uniqid from 'uniqid';

const registeredOnStart = new Map();
const registeredOnCompleted = new Map();

if (is.renderer) {
  ipcRenderer.removeAllListeners('download-started');
  ipcRenderer.on('download-started', (event, { id }) => {
    if (registeredOnStart.has(id)) {
      const onStart = registeredOnStart.get(id);
      if (typeof onStart === 'function') {
        onStart();
      }
    }
  });
  ipcRenderer.removeAllListeners('download-completed');
  ipcRenderer.on('download-completed', (event, { id }) => {
    if (registeredOnCompleted.has(id)) {
      const onCompleted = registeredOnCompleted.get(id);
      if (typeof onCompleted === 'function') {
        onCompleted();
      }
      registeredOnCompleted.delete(id);
    }
    if (registeredOnStart.has(id)) {
      registeredOnStart.delete(id);
    }
  });
}

export default {
  downloadUrl(
    url: string,
    filename: string,
    onStart?: () => void,
    onCompleted?: () => void
  ) {
    const id = uniqid();
    if (onStart) registeredOnStart.set(id, onStart);
    if (onCompleted) registeredOnCompleted.set(id, onCompleted);
    ipcRenderer.send('download-file', { id, url, filename });
  }
};
