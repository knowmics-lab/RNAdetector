// @flow
/* eslint-disable no-restricted-syntax */
import fs from 'fs-extra';
import { is } from 'electron-util';
import Client from 'dockerode';
import Utils from './utils';
import Settings from './settings';
import type { ConfigObjectType } from '../types/settings';
import type { Package } from '../types/local';

export const DOCKER_IMAGE_NAME = 'alaimos/rnadetector:v0.0.1';

type DockerPullEvent = {
  status: string,
  id?: string,
  progress?: string
};

export class DockerPullStatus {
  idMap: Map<string, number> = new Map<string, number>();

  outputArray: string[] = [];

  pushEvent(event: DockerPullEvent) {
    if (event.id) {
      const { id, status } = event;
      let mappedId;
      if (this.idMap.has(id)) {
        mappedId = this.idMap.get(id);
      } else {
        mappedId = this.outputArray.length;
        this.idMap.set(id, mappedId);
        this.outputArray.push('');
      }
      if (mappedId) {
        const progress = event.progress ? ` ${event.progress}` : '';
        this.outputArray[mappedId] = `${id}: ${status}${progress}`;
      }
    } else {
      this.outputArray.push(event.status);
    }
  }

  toString() {
    return this.outputArray.join('\n');
  }
}

export class DockerManager {
  config: ConfigObjectType;

  client: Client;

  container: Client.Container;

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
    stream: http$IncomingMessage<>,
    onStdout: ?(Buffer) => void,
    onStderr: ?(Buffer) => void,
    onEnd: ?() => void
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

  static async demuxStream(
    stream: http$IncomingMessage<>
  ): Promise<[string, string]> {
    return new Promise(resolve => {
      let stdout = Buffer.from('');
      let stderr = Buffer.from('');
      DockerManager.liveDemuxStream(
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
    try {
      const inspect = await this.getContainer().inspect();
      return inspect.State.Status;
    } catch (e) {
      if (e.statusCode && e.statusCode === 404) {
        return 'not found';
      }
      throw e;
    }
  }

  getContainer() {
    if (!this.container) {
      this.container = this.client.getContainer(this.config.containerName);
      /* const containers = await this.client.listContainers({
        all: true
      });
      const found = containers
        .filter(c => c.Image === DOCKER_IMAGE_NAME)
        .filter(c => c.Names.includes(`/${this.config.containerName}`));
      if (found.length === 0) {
        return null;
      }
      this.container = this.client.getContainer(found[0].Id); */
    }
    return this.container;
  }

  async createContainer() {
    const status = await this.checkContainerStatus();
    if (status === 'not found') {
      await this.cleanupBootedFile();
      this.container = null;
      // noinspection ES6MissingAwait
      this.client.createContainer({
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
      return new Promise((resolve, reject) => {
        const timer = setInterval(async () => {
          const currentStatus = await this.checkContainerStatus();
          if (currentStatus !== 'not found') {
            clearInterval(timer);
            if (currentStatus !== 'created') {
              reject(
                new Error(
                  `Unable to create the container ${this.config.containerName}. Create it manually`
                )
              );
            }
            try {
              await this.startContainer();
              resolve();
            } catch (e) {
              reject(e);
            }
          }
        }, 500);
      });
    }
  }

  async startContainer() {
    const status = await this.checkContainerStatus();
    if (status === 'not found') return this.createContainer();
    if (status === 'exited' || status === 'created') {
      await this.cleanupBootedFile();
      const container = this.getContainer();
      // noinspection ES6MissingAwait
      container.start();
      return new Promise((resolve, reject) => {
        const timer = setInterval(async () => {
          const currentStatus = await this.checkContainerStatus();
          if (currentStatus !== 'exited' && currentStatus !== 'created') {
            clearInterval(timer);
            if (currentStatus !== 'running') {
              reject(
                new Error(
                  `Unable to start the container ${this.config.containerName}. Start it manually`
                )
              );
            }
            try {
              await this.waitContainerBooted();
              resolve();
            } catch (e) {
              reject(e);
            }
          }
        }, 500);
      });
    }
  }

  async stopContainer() {
    const status = await this.checkContainerStatus();
    if (status === 'running') {
      const container = this.getContainer();
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
    const container = this.getContainer();
    if (container) {
      await container.remove();
    } else {
      throw new Error('Unable to find container');
    }
  }

  async execDockerCommand(Cmd: string[]): Promise<*> {
    const status = await this.checkContainerStatus();
    if (status === 'running') {
      const container = this.getContainer();
      if (!container) throw new Error('Unable to get container instance');
      const exec = await container.exec({
        Cmd,
        AttachStdout: true
      });
      const stream = await exec.start();
      const [stdout] = await DockerManager.demuxStream(stream);
      console.log(stdout);
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

  async execDockerCommandLive(
    Cmd: string[],
    outputCallback: string => void,
    errCallback: ?(string) => void = null,
    exitCallback: ?(number) => void = null
  ): Promise<void> {
    const status = await this.checkContainerStatus();
    if (status === 'running') {
      const container = this.getContainer();
      if (!container) throw new Error('Unable to get container instance');
      const exec = await container.exec({
        Cmd,
        AttachStdout: true,
        AttachStderr: !!errCallback
      });
      const stream = await exec.start();
      const onStderr = errCallback ? buf => errCallback(buf.toString()) : null;
      const onExit = exitCallback
        ? () => {
            exec.inspect((err, data) => {
              if (err) throw new Error(err);
              exitCallback(data.ExitCode);
            });
          }
        : null;
      return DockerManager.liveDemuxStream(
        stream,
        buf => outputCallback(buf.toString()),
        onStderr,
        onExit
      );
    }
    throw new Error('Unable to exec command. Container is not running');
  }

  installPackage(
    name: string,
    outputCallback: string => void,
    errorCallback: (*) => void,
    exitCallback: number => void
  ): void {
    this.execDockerCommandLive(
      ['php', '/rnadetector/ws/artisan', 'packages:install', name],
      outputCallback,
      outputCallback,
      exitCallback
    ).catch(e => errorCallback(e));
  }

  async pullImage(outputCallback: ?(DockerPullStatus) => void) {
    return new Promise((resolve, reject) => {
      this.client.pull(DOCKER_IMAGE_NAME, (e, stream) => {
        if (e) {
          reject(e);
        } else {
          const status = new DockerPullStatus();
          const onFinished = err => {
            if (err) reject(err);
            else resolve(status);
          };
          const onProgress = event => {
            if (outputCallback) {
              status.pushEvent(event);
              outputCallback(status);
            }
          };
          this.client.modem.followProgress(stream, onFinished, onProgress);
        }
      });
    });
  }

  async hasImage() {
    const images = await this.client.listImages();
    return (
      images.filter(r => r.RepoTags.includes(DOCKER_IMAGE_NAME)).length > 0
    );
  }
}

let instance = null;

export const getInstance = () => {
  if (!instance) {
    instance = new DockerManager(Settings.getConfig());
    if (is.renderer) window.docker = instance;
  }
  return instance;
};

export const resetInstance = () => {
  instance = null;
};

export default getInstance();
