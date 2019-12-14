import Store from 'electron-store';
import configSchema from '../constants/config-schema.json';

export type ConfigObjectType = {
  +local: boolean,
  +webserviceUrl: string,
  +jobsPath: string,
  +containerName: string,
  +apiKey: string
};

export default {
  configStore: new Store({ schema: configSchema }),
  getConfig() {
    return {
      local: this.configStore.get('local'),
      webserviceUrl: this.configStore.get('webserviceUrl'),
      jobsPath: this.configStore.get('jobsPath'),
      containerName: this.configStore.get('containerName'),
      apiKey: this.configStore.get('apiKey')
    };
  },
  saveConfig(config: ConfigObjectType) {
    return new Promise((resolve, reject) => {
      try {
        this.configStore.set(config);
        resolve(config);
      } catch (e) {
        reject(e.message);
      }
    });
  }
};
