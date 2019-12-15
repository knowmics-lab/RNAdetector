// @flow
import fs from 'fs';
import path from 'path';

let watcher = null;

export default {
  async waitExists(filePath: string, timeout: number = 0): Promise<*> {
    return new Promise((resolve, reject) => {
      let timer = null;
      const closeWatcher = () => {
        if (watcher !== null) {
          watcher.close();
          watcher = null;
        }
      };
      const closeTimeout = () => {
        if (timer !== null) {
          clearTimeout(timer);
        }
      };
      if (timeout > 0) {
        timer = setTimeout(() => {
          closeWatcher();
          reject(new Error(`Unable to find ${filePath}. Operation timed out`));
        }, timeout);
      }

      fs.access(filePath, fs.constants.R_OK, err => {
        if (!err) {
          closeTimeout();
          closeWatcher();
          resolve();
        }
      });

      const dir = path.dirname(filePath);
      const basename = path.basename(filePath);
      watcher = fs.watch(dir, (eventType, filename) => {
        if (eventType === 'rename' && filename === basename) {
          closeTimeout();
          closeWatcher();
          resolve();
        }
      });
    });
  }
};
