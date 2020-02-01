// @flow
import fs from 'fs';
import path from 'path';
import os from 'os';
import type { FileFilter } from '../types/common';
import type { AnalysisFileTypes } from '../types/analysis';

let watcher = null;

export default {
  cpuCount() {
    return os.cpus().length;
  },
  supportedAnalysisFileTypes() {
    return {
      fastq: 'FASTQ',
      bam: 'BAM',
      sam: 'SAM'
    };
  },
  analysisFileExtensions(type: AnalysisFileTypes): FileFilter[] {
    switch (type) {
      case 'fastq':
        return [
          { name: 'FASTQ files', extensions: ['fq', 'fastq', 'gz', 'bz2'] }
        ];
      case 'bam':
        return [{ name: 'BAM files', extensions: ['bam'] }];
      case 'sam':
        return [{ name: 'SAM files', extensions: ['sam'] }];
      default:
        return [];
    }
  },
  filterByKey(raw: {}, callback: string => boolean) {
    return Object.keys(raw)
      .filter(callback)
      .reduce((obj, key) => {
        return {
          ...obj,
          [key]: raw[key]
        };
      }, {});
  },
  toArray(list: *) {
    return Array.prototype.slice.call(list || [], 0);
  },
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
  },
  dashToWordString(s: string) {
    return s.replace(/[_\\-]([a-z0-9])/g, g => ` ${g[1].toUpperCase()}`);
  },
  capitalize(s: string) {
    return s.charAt(0).toUpperCase() + s.slice(1);
  }
};
