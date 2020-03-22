// @flow
/* eslint-disable no-restricted-syntax */
import process from 'child_process';
import fs from 'fs-extra';
import util from 'util';
import { is } from 'electron-util';
import Client from 'dockerode';
import Utils from './utils';
// eslint-disable-next-line import/no-cycle
import Settings from './settings';
import type { ConfigObjectType } from '../types/settings';
import type { Package } from '../types/local';

const execFile = util.promisify(process.execFile);

export const DOCKER_IMAGE_NAME = 'alaimos/rnadetector:v0.0.1';
export const VALID_DOCKER_STATE = [
  'created',
  'running',
  'paused',
  'restarting',
  'removing',
  'exited',
  'dead'
];

class Docker {
  config: ConfigObjectType;

  client;

  container;

  constructor(config: ConfigObjectType) {
    this.config = config;
    let socket;
    if (config.socketPath) {
      socket = config.socketPath;
    } else if (is.windows) {
      socket = '//./pipe/docker_engine';
    } else {
      socket = '/var/run/docker.sock';
    }
    this.client = new Client({
      socketPath: socket
    });
  }

  static liveDemuxStream(
    stream,
    onStdout: ?(Buffer) => void,
    onStderr: ?(Buffer) => void,
    onEnd: () => void
  ): void {
    let nextDataType = null;
    let nextDataLength = -1;
    let buffer = Buffer.from('');

    const bufferSlice = (end: number) => {
      const out = buffer.slice(0, end);
      buffer = Buffer.from(buffer.slice(end, buffer.length));
      return out;
    };
    const processData = data => {
      if (data) {
        buffer = Buffer.concat([buffer, data]);
      }
      if (nextDataType) {
        if (buffer.length >= nextDataLength) {
          const content = bufferSlice(nextDataLength);
          if (onStdout && nextDataType === 1) {
            onStdout(Buffer.from(content));
          } else if (onStderr && nextDataType !== 1) {
            onStderr(Buffer.from(content));
          }
          nextDataType = null;
          processData();
        }
      } else if (buffer.length >= 8) {
        const header = bufferSlice(8);
        nextDataType = header.readUInt8(0);
        nextDataLength = header.readUInt32BE(4);
        processData();
      }
    };

    stream.on('data', processData);
    if (onEnd) {
      stream.on('end', onEnd);
    }
  }

  static async demuxStream(stream): Promise<[string, string]> {
    return new Promise(resolve => {
      let stdout = Buffer.from('');
      let stderr = Buffer.from('');
      Docker.liveDemuxStream(
        stream,
        content => {
          stdout = Buffer.concat([stdout, content]);
        },
        content => {
          stderr = Buffer.concat([stderr, content]);
        },
        () => resolve([stdout.toString(), stderr.toString()])
      );
    });
  }

  getBootedFile(): string {
    return `${this.config.dataPath}/booted`;
  }

  getDbReadyFile(): string {
    return `${this.config.dataPath}/database/ready`;
  }

  async waitContainerBooted() {
    await Utils.waitExists(this.getDbReadyFile());
    await Utils.waitExists(this.getBootedFile());
  }

  async cleanupBootedFile() {
    await fs.remove(this.getBootedFile());
  }

  async checkContainerStatus() {
    const containers = await this.client.listContainers({
      all: true
    });
    const found = containers
      .filter(c => c.Image === DOCKER_IMAGE_NAME)
      .filter(c => c.Names.includes(`/${this.config.containerName}`));
    if (found.length === 0) {
      return 'not found';
    }
    if (VALID_DOCKER_STATE.includes(found[0].State)) {
      return found[0].State;
    }
    throw new Error('unknown state');
  }

  async getContainer() {
    if (!this.container) {
      const containers = await this.client.listContainers({
        all: true
      });
      const found = containers
        .filter(c => c.Image === DOCKER_IMAGE_NAME)
        .filter(c => c.Names.includes(`/${this.config.containerName}`));
      if (found.length === 0) {
        return null;
      }
      this.container = this.client.getContainer(found[0].Id);
    }
    return this.container;
  }

  async createContainer() {
    const status = await this.checkContainerStatus();
    if (status === 'not found') {
      await this.cleanupBootedFile();
      const container = await this.client.createContainer({
        Image: DOCKER_IMAGE_NAME,
        name: this.config.containerName,
        ExposedPorts: {
          '80/tcp': {}
        },
        Volumes: {
          '/rnadetector/ws/storage/app/': {}
        },
        HostConfig: {
          PortBindings: {
            '80/tcp': [
              {
                HostPort: `${this.config.apiPort}`
              }
            ]
          },
          Binds: [`${this.config.dataPath}:/rnadetector/ws/storage/app/`]
        }
      });
      await container.start();
      await this.waitContainerBooted();
    }
  }

  async startContainer() {
    const status = await this.checkContainerStatus();
    if (status === 'not found') return this.createContainer();
    if (status === 'exited') {
      await this.cleanupBootedFile();
      const container = await this.getContainer();
      if (container) {
        await container.start();
        if ((await this.checkContainerStatus()) !== 'running') {
          throw new Error(
            `Unable to start the container ${this.config.containerName}. Start it manually`
          );
        }
        await this.waitContainerBooted();
      } else {
        throw new Error('Unable to find container');
      }
    }
  }

  async stopContainer() {
    const status = await this.checkContainerStatus();
    if (status === 'running') {
      const container = await this.getContainer();
      if (container) {
        await container.stop();
        await this.cleanupBootedFile();
      } else {
        throw new Error('Unable to find container');
      }
    }
  }

  async removeContainer() {
    const status = await this.checkContainerStatus();
    if (status === 'running') {
      await this.stopContainer();
    }
    const container = await this.getContainer();
    if (container) {
      await container.remove();
    } else {
      throw new Error('Unable to find container');
    }
  }

  async execDockerCommand(Cmd: string[]): Promise<*> {
    const status = await this.checkContainerStatus();
    if (status === 'running') {
      const container = await this.getContainer();
      if (!container) throw new Error('Unable to get container instance');
      const exec = await container.exec({
        Cmd,
        AttachStdout: true
      });
      const stream = await exec.start();
      const [stdout] = await Docker.demuxStream(stream);
      return JSON.parse(stdout);
    }
    throw new Error('Unable to exec command. Container is not running');
  }

  async generateAuthToken(): Promise<string> {
    const result = await this.execDockerCommand(['/genkey.sh', '--json']);
    if (!result.error) {
      return result.data;
    }
    throw new Error(result.message);
  }

  async listPackages(): Promise<Package[]> {
    const result = await this.execDockerCommand([
      'php',
      '/rnadetector/ws/artisan',
      'packages:list'
    ]);
    if (!result.error) {
      return result.packages;
    }
    throw new Error(result.message);
  }
}

export default {
  getBootedFile(config: ConfigObjectType = Settings.getConfig()): string {
    return `${config.dataPath}/booted`;
  },
  getDbReadyFile(config: ConfigObjectType = Settings.getConfig()): string {
    return `${config.dataPath}/database/ready`;
  },
  async waitContainerBooted(config: ConfigObjectType = Settings.getConfig()) {
    await Utils.waitExists(this.getDbReadyFile(config));
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
    if (is.renderer && window && !window.docker)
      window.docker = new Docker(Settings.getConfig());
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
  },
  async listPackages(
    config: ConfigObjectType = Settings.getConfig()
  ): Promise<Package[]> {
    const result = await this.execDockerCommand(
      ['php', '/rnadetector/ws/artisan', 'packages:list'],
      config
    );
    if (!result.error) {
      return result.packages;
    }
    throw new Error(result.message);
  },
  async execDockerCommandLive(
    command: string[],
    outputCallback: string => void,
    errCallback: ?(string) => void = null,
    exitCallback: ?(number) => void = null,
    config: ConfigObjectType = Settings.getConfig()
  ): Promise<void> {
    const status = await this.checkContainerStatus(config);
    if (status === 'running') {
      const child = process.execFile(config.dockerExecutablePath, [
        'exec',
        config.containerName,
        ...command
      ]);
      child.stdout.on('data', data => outputCallback(data.toString()));
      if (errCallback)
        child.stderr.on('data', data => errCallback(data.toString()));
      if (exitCallback) child.on('exit', code => exitCallback(code));
    } else {
      throw new Error('Unable to exec command. Container is not running');
    }
  },
  installPackage(
    name: string,
    outputCallback: string => void,
    errorCallback: (*) => void,
    exitCallback: number => void,
    config: ConfigObjectType = Settings.getConfig()
  ): void {
    this.execDockerCommandLive(
      ['php', '/rnadetector/ws/artisan', 'packages:install', name],
      outputCallback,
      null,
      exitCallback,
      config
    ).catch(e => errorCallback(e));
  }
};
