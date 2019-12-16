// @flow
/* eslint-disable no-restricted-syntax */
import process from 'child_process';
import fs from 'fs-extra';
import util from 'util';
import Utils from './utils';
// eslint-disable-next-line import/no-cycle
import Settings from './settings';
import type { ConfigObjectType } from './settings';

const execFile = util.promisify(process.execFile);

export const DOCKER_IMAGE_NAME = 'alaimos/ubuntu-private:RNAdetector.v1.0';

export default {
  getBootedFile(config: ConfigObjectType = Settings.getConfig()): string {
    return `${config.dataPath}/booted`;
  },
  async waitContainerBooted(config: ConfigObjectType = Settings.getConfig()) {
    await Utils.waitExists(this.getBootedFile(config));
  },
  async cleanupBootedFile(config: ConfigObjectType = Settings.getConfig()) {
    await fs.remove(this.getBootedFile(config));
  },
  checkDockerProcess(
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
    const status = await this.checkContainerStatus(config);
    if (status === 'not found') {
      await this.cleanupBootedFile();
      await execFile(config.dockerExecutablePath, [
        'run',
        '-d',
        '-p',
        `${config.apiPort}:80`,
        '-v',
        `${config.dataPath}:/rnadetector/ws/storage/app/`,
        `--name=${config.containerName}`,
        DOCKER_IMAGE_NAME
      ]);
      if ((await this.checkContainerStatus(config)) !== 'running') {
        throw new Error(
          `Unable to create the container ${config.containerName}. Create it manually`
        );
      }
      await this.waitContainerBooted(config);
    }
  },
  async startContainer(config: ConfigObjectType = Settings.getConfig()) {
    const status = await this.checkContainerStatus(config);
    if (status === 'not found') return this.createContainer(config);
    if (status === 'stopped') {
      await this.cleanupBootedFile();
      await execFile(config.dockerExecutablePath, [
        'start',
        config.containerName
      ]);
      if ((await this.checkContainerStatus(config)) !== 'running') {
        throw new Error(
          `Unable to start the container ${config.containerName}. Start it manually`
        );
      }
      await this.waitContainerBooted(config);
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
      await this.cleanupBootedFile();
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
  },
  async execDockerCommand(
    command: string[],
    config: ConfigObjectType = Settings.getConfig()
  ): Promise<*> {
    const status = await this.checkContainerStatus(config);
    if (status === 'running') {
      const { stdout } = await execFile(config.dockerExecutablePath, [
        'exec',
        config.containerName,
        ...command
      ]);
      return JSON.parse(stdout);
    }
    throw new Error('Unable to exec command. Container is not running');
  },
  async generateAuthToken(
    config: ConfigObjectType = Settings.getConfig()
  ): Promise<string> {
    const result = await this.execDockerCommand(
      ['/genkey.sh', '--json'],
      config
    );
    if (!result.error) {
      return result.data;
    }
    throw new Error(result.message);
  }
};
