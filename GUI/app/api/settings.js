// @flow
import Store from 'electron-store';
import { api } from 'electron-util';
import fs from 'fs-extra';
import axios from 'axios';
import configSchema from '../constants/config-schema.json';
// eslint-disable-next-line import/no-cycle
import Docker from './docker';

export type ConfigObjectType = {
  +local: boolean,
  +webserviceUrl: string,
  +jobsPath: string,
  +containerName: string,
  +apiKey: string,
  +dockerExecutablePath: string
};

export default {
  configStore: new Store({ schema: configSchema }),
  getConfig(): ConfigObjectType {
    return {
      local: this.configStore.get('local'),
      webserviceUrl: this.configStore.get('webserviceUrl'),
      jobsPath: this.configStore.get(
        'jobsPath',
        `${api.app.getPath('home')}/.RNADetector`
      ),
      containerName: this.configStore.get('containerName'),
      apiKey: this.configStore.get('apiKey'),
      dockerExecutablePath: this.configStore.get('dockerExecutablePath')
    };
  },
  saveConfig(config: ConfigObjectType): Promise<ConfigObjectType> {
    // eslint-disable-next-line no-async-promise-executor
    return new Promise(async (resolve, reject) => {
      try {
        const oldConfig = this.getConfig();
        const newConfig = {
          ...oldConfig,
          ...config
        };
        if (oldConfig.containerName !== newConfig.containerName) {
          Docker.removeContainer(oldConfig);
        }
        await this.checkConfig(newConfig);
        this.configStore.set(newConfig);
        resolve(newConfig);
      } catch (e) {
        reject(e.message);
      }
    });
  },
  async checkUrl(config: ConfigObjectType) {
    const {
      data: { data }
    } = await axios.get(`${config.webserviceUrl}/api/ping`);
    if (data !== 'pong') throw new Error('Invalid webservice URL');
  },
  async checkToken(config: ConfigObjectType) {
    let data = null;
    try {
      const response = await axios.get(`${config.webserviceUrl}/api/auth-ping`);
      data = response.data.data;
    } catch (e) {
      throw new Error(`Invalid authentication token - ${e.message}`);
    }
    if (data !== 'pong') throw new Error('Invalid authentication token');
  },
  async checkConfig(config: ConfigObjectType = this.getConfig()) {
    if (config.local) {
      if (!(await fs.pathExists(config.jobsPath))) {
        await fs.ensureDir(config.jobsPath, 0o755);
      }
    }
    await Docker.checkDockerProcess(config);
    const status = await Docker.checkContainerStatus(config);
    if (status !== 'running') {
      await Docker.startContainer(config);
    }
    await this.checkUrl(config);
    // await this.checkToken(config);
  }
};
