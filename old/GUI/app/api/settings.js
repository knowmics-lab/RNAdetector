// @flow
import Store from 'electron-store';
import { api } from 'electron-util';
import findFreePort from 'find-free-port';
import configSchema from '../constants/config-schema.json';
import type { ConfigObjectType } from '../types/settings';
import type { AxiosHeaders, SimpleMapType } from '../types/common';

export default {
  configStore: new Store({ schema: configSchema }),
  getConfig(): ConfigObjectType {
    return {
      configured: this.configStore.get('configured'),
      local: this.configStore.get('local'),
      apiProtocol: this.configStore.get('apiProtocol'),
      apiHostname: this.configStore.get('apiHostname'),
      apiPort: this.configStore.get('apiPort'),
      apiPath: this.configStore.get('apiPath'),
      publicPath: this.configStore.get('publicPath'),
      dataPath: this.configStore.get(
        'dataPath',
        `${api.app.getPath('home')}/.RNADetector`
      ),
      socketPath: this.configStore.get('socketPath'),
      containerName: this.configStore.get('containerName'),
      apiKey: this.configStore.get('apiKey'),
      autoStopDockerOnClose: this.configStore.get('autoStopDockerOnClose')
    };
  },
  autoStopDockerOnClose() {
    return this.configStore.get('autoStopDockerOnClose');
  },
  setAutoStopDockerOnClose() {
    const newConfig = {
      ...this.getConfig(),
      autoStopDockerOnClose: true
    };
    this.configStore.set(newConfig);
  },
  isConfigured() {
    return this.configStore.get('configured');
  },
  isLocal() {
    return this.configStore.get('local');
  },
  getApiUrl(): string {
    const config = this.getConfig();
    const path = config.apiPath.replace(/^\/|\/$/gm, '');
    return `${config.apiProtocol}://${config.apiHostname}:${config.apiPort}/${path}/`;
  },
  getJBrowseUrl(configUri: string): string {
    const config = this.getConfig();
    return `${config.apiProtocol}://${config.apiHostname}:${
      config.apiPort
    }/jbrowse2/index.html?config=${encodeURIComponent(configUri)}`;
  },
  getPublicUri(p: string = ''): string {
    const config = this.getConfig();
    const path = config.publicPath.replace(/^\/|\/$/gm, '');
    return `/${path}/${p ? p.replace(/^\//, '') : ''}`;
  },
  getPublicUrl(p: string = ''): string {
    const config = this.getConfig();
    const path = config.publicPath.replace(/^\/|\/$/gm, '');
    return `${config.apiProtocol}://${config.apiHostname}:${
      config.apiPort
    }/${path}/${p ? p.replace(/^\//, '') : ''}`;
  },
  getLocalPath(p: string = '') {
    const config = this.getConfig();
    return `${config.dataPath}/${p ? p.replace(/^\//, '') : ''}`;
  },
  getAuthHeaders(): SimpleMapType<string> {
    const config = this.getConfig();
    return {
      Authorization: `Bearer ${config.apiKey}`
    };
  },
  getAxiosHeaders(): AxiosHeaders {
    const config = this.getConfig();
    return {
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
        Authorization: `Bearer ${config.apiKey}`
      }
    };
  },
  saveConfig(config: ConfigObjectType) {
    return this.configStore.set(config);
  },
  async findFreePort(start: number): Promise<number> {
    const [port] = await findFreePort(start);
    return port;
  }
};
