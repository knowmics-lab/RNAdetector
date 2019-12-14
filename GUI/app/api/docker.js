// @flow
/* eslint-disable no-restricted-syntax */
import process from 'child_process';
import { is } from 'electron-util';
import util from 'util';
// eslint-disable-next-line import/no-cycle
import Settings from './settings';
import type { ConfigObjectType } from './settings';

const execFile = util.promisify(process.execFile);

export default {
  async checkDockerProcess(
    config: ConfigObjectType = Settings.getConfig()
  ): Promise<void> {
    return new Promise((resolve, reject) => {
      const commonReject = () => reject(new Error('Invalid docker executable'));
      process
        .spawn(config.dockerExecutablePath, ['version'])
        .on('error', commonReject)
        .on('close', code => {
          if (code !== 0) commonReject();
          resolve();
        });
    });
  },
  async checkContainerStatus(config: ConfigObjectType = Settings.getConfig()) {
    const { stdout } = await execFile(config.dockerExecutablePath, [
      'ps',
      '-a',
      '-f',
      `name=${config.containerName}`,
      '--format',
      "'{{.Names}}\\t{{.Status}}'"
    ]);
    if (!stdout) return 'not found';
    for (const r of stdout.replace(/^'|'$/gm, '').split('\n')) {
      const [name, status] = r.split('\t');
      if (name === config.containerName) {
        if (status.toLowerCase().includes('exited')) {
          return 'stopped';
        }
        if (status.toLowerCase().includes('up')) {
          return 'running';
        }
        throw new Error(`Unable to parse status ${status}`);
      }
    }
    return 'not found';
  },
  async createContainer(config: ConfigObjectType = Settings.getConfig()) {
    // docker run -d -p 8888:80 -v /home/alaimos/rnadet/:/rnadetector/ws/storage/app/ --name=RNAdetector alaimos/ubuntu-private:RNAdetector.v1.0
  },
  async startContainer(config: ConfigObjectType = Settings.getConfig()) {
    const status = await this.checkContainerStatus(config);
    if (status === 'not found') this.createContainer(config);
    if (status === 'stopped') {
      await execFile(config.dockerExecutablePath, [
        'start',
        config.containerName
      ]);
      if ((await this.checkContainerStatus(config)) !== 'running') {
        throw new Error(
          `Unable to start the container ${config.containerName}. Start it manually`
        );
      }
    }
  },
  async stopContainer(config: ConfigObjectType = Settings.getConfig()) {
    const status = await this.checkContainerStatus(config);
    if (status === 'running') {
      await execFile(config.dockerExecutablePath, [
        'stop',
        config.containerName
      ]);
      if ((await this.checkContainerStatus(config)) !== 'stopped') {
        throw new Error(
          `Unable to stop the container ${config.containerName}. Stop it manually`
        );
      }
    }
  },
  async removeContainer(config: ConfigObjectType = Settings.getConfig()) {
    const status = await this.checkContainerStatus(config);
    if (status === 'running') this.stopContainer(config);
    await execFile(config.dockerExecutablePath, ['rm', config.containerName]);
    if ((await this.checkContainerStatus(config)) !== 'not found') {
      throw new Error(
        `Unable to remove the container ${config.containerName}. Remove it manually`
      );
    }
  }
};
