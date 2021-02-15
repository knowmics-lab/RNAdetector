// @flow
/* eslint-disable no-restricted-syntax,max-classes-per-file */
import fs from 'fs-extra';
import { is } from 'electron-util';
import Client from 'dockerode';
import Utils from './utils';
import Settings from './settings';
import { DOCKER_IMAGE_NAME } from '../constants/system.json';
import type { ConfigObjectType } from '../types/settings';
import type { Package } from '../types/local';
import TimeoutError from '../errors/TimeoutError';

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

  isUpToDate() {
    return (
      this.outputArray.filter(s => s.includes('Image is up to date')).length > 0
    );
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
    this.initClient();
  }

  initClient() {
    let socket;
    if (this.config.socketPath) {
      socket = this.config.socketPath;
    } else if (is.windows) {
      socket = '//./pipe/docker_engine';
    } else {
      socket = '/var/run/docker.sock';
    }
    this.client = new Client({
      socketPath: socket
    });
    this.container = undefined;
  }

  static liveDemuxStream(
    stream: http$IncomingMessage<>,
    onStdout: ?(Buffer) => void,
    onStderr: ?(Buffer) => void,
    onEnd: ?() => void,
    checkRunning: ?() => Promise<boolean>,
    timeoutRunning: number = 30000
  ): void {
    let nextDataType = null;
    let nextDataLength = -1;
    let buffer = Buffer.from('');
    let ended = false;

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

    stream.on('data', processData).on('end', () => {
      if (!ended && onEnd) {
        onEnd();
        ended = true;
      }
    });
    if (checkRunning) {
      const fnRunning = async () => {
        if (ended) return;
        if (await checkRunning()) {
          setTimeout(fnRunning, timeoutRunning);
        } else if (!ended && onEnd) {
          onEnd();
          ended = true;
        }
      };
      setTimeout(fnRunning, timeoutRunning);
    }
  }

  static async demuxStream(
    stream: http$IncomingMessage<>,
    checkRunning: ?() => Promise<boolean>,
    timeoutRunning: number = 30000
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
        () => resolve([stdout.toString(), stderr.toString()]),
        checkRunning,
        timeoutRunning
      );
    });
  }

  getBootedFile(): string {
    return `${this.config.dataPath}/booted`;
  }

  getDbDirectory(): string {
    return `${this.config.dataPath}/database/`;
  }

  getDbReadyFile(): string {
    return `${this.getDbDirectory()}/ready`;
  }

  async waitContainerBooted(timeout: number = 0) {
    await Utils.waitExists(this.getDbDirectory(), timeout);
    await Utils.waitExists(this.getDbReadyFile(), timeout);
    await Utils.waitExists(this.getBootedFile(), timeout);
  }

  async cleanupBootedFile() {
    await fs.remove(this.getBootedFile());
  }

  async checkContainerStatus(): Promise<string> {
    let inspect = null;
    while (inspect === null) {
      try {
        // eslint-disable-next-line no-await-in-loop
        inspect = await Utils.promiseTimeout(
          this.getContainer().inspect(),
          500
        );
      } catch (e) {
        if (!(e instanceof TimeoutError)) {
          if (e.statusCode && e.statusCode === 404) {
            return 'not found';
          }
          throw e;
        }
      }
    }
    return inspect.State.Status;
  }

  async isRunning() {
    const status = await this.checkContainerStatus();
    return status === 'running';
  }

  getContainer() {
    if (!this.container) {
      this.container = this.client.getContainer(this.config.containerName);
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

  async checkForUpdates(
    showMessage: (string, boolean) => void,
    displayLog: ?(string) => void,
    timeout: number = 120000,
    maxTries: number = 3
  ): Promise<boolean> {
    if (this.config.local) {
      showMessage('Checking internet connection...', false);
      if (!(await this.isRunning()) && (await Utils.isOnline())) {
        showMessage('Checking for container updates...', false);
        try {
          let status = null;
          let timer;
          const displayStatus = displayLog
            ? ns => {
                status = ns;
              }
            : undefined;
          if (displayLog) {
            timer = setInterval(() => {
              if (status) displayLog(status.toString());
              status = null;
            }, 500);
          }
          const res = await this.pullImage(displayStatus);
          if (displayLog && timer) {
            clearInterval(timer);
          }
          if (!res.isUpToDate()) {
            await Utils.retryFunction(
              async (t: number) => {
                const first = t === 0;
                showMessage(
                  `Update found...removing old container...${
                    first ? '' : `Attempt ${t + 1} of ${maxTries}...`
                  }`,
                  !first
                );
                await this.removeContainer();
              },
              timeout,
              maxTries
            );
            return true;
          }
        } catch (e) {
          if (e instanceof TimeoutError) {
            throw new Error(
              `Unable to update the container ${this.config.containerName}. Update it manually`
            );
          } else {
            throw e;
          }
        }
      }
    }
    return false;
  }

  async runUpdateScript(
    showMessage: (string, boolean) => void,
    displayLog: ?(string) => void
  ): Promise<void> {
    if (await this.isUpdateScriptNeeded()) {
      showMessage('Running update script...', false);
      let log = '';
      let timer;
      const displayStatus = l => {
        log += l;
      };
      if (displayLog) {
        timer = setInterval(() => {
          if (log) displayLog(log);
        }, 500);
      }
      return new Promise((resolve, reject) => {
        this.execDockerCommandLive(
          ['/update_run.sh'],
          displayStatus,
          displayStatus,
          code => {
            if (code !== 0) {
              showMessage(`An unknown error occurred! Code: ${code}.`, true);
            }
            if (displayLog && timer) {
              clearInterval(timer);
            }
            resolve();
          }
        ).catch(reject);
      });
    }
  }

  async startupSequence(
    showMessage: (string, boolean) => void,
    displayLog: ?(string) => void,
    timeout: number = 120000,
    maxTries: number = 3
  ) {
    const updated = await this.checkForUpdates(
      showMessage,
      displayLog,
      timeout,
      maxTries
    );
    try {
      await Utils.retryFunction(
        async (t: number) => {
          const first = t === 0;
          const odd = t % 2 > 0;
          showMessage(
            first
              ? 'Starting docker container...'
              : `Container is not starting...Attempt ${t +
                  1} of ${maxTries}...`,
            !first
          );
          if (odd) {
            await this.removeContainer();
          }
          return this.startContainer();
        },
        timeout,
        maxTries
      );
      if (updated) {
        await this.runUpdateScript(showMessage, displayLog);
      }
    } catch (e) {
      if (e instanceof TimeoutError) {
        throw new Error(
          `Unable to start the container ${this.config.containerName}. Start it manually`
        );
      } else {
        throw e;
      }
    }
  }

  async startContainer() {
    const status = await this.checkContainerStatus();
    if (status === 'not found') return this.createContainer();
    if (status === 'exited' || status === 'created') {
      await this.cleanupBootedFile();
      const container = this.getContainer();
      container.start();
      return new Promise((resolve, reject) => {
        let timer;
        const fnTimer = async () => {
          const currentStatus = await this.checkContainerStatus();
          if (currentStatus === 'running') {
            timer = clearInterval(timer);
            try {
              await this.waitContainerBooted();
              return resolve();
            } catch (e) {
              return reject(e);
            }
          }
        };
        timer = setInterval(fnTimer, 500);
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
    if (status !== 'not found') {
      const container = this.getContainer();
      if (container) {
        await container.remove();
      } else {
        throw new Error('Unable to find container');
      }
    }
  }

  async execDockerCommand(
    Cmd: string[],
    timeoutRunning: number = 30000,
    parse: boolean = true
  ): Promise<*> {
    const status = await this.checkContainerStatus();
    if (status === 'running') {
      const container = this.getContainer();
      if (!container) throw new Error('Unable to get container instance');
      const exec = await container.exec({
        Cmd,
        AttachStdout: true
      });
      const stream = await exec.start();
      const [stdout] = await DockerManager.demuxStream(
        stream,
        () => {
          return new Promise(resolve => {
            exec.inspect((e, d) => {
              resolve(d && d.Running);
            });
          });
        },
        timeoutRunning
      );
      if (parse) {
        return JSON.parse(stdout);
      }
      return stdout;
    }
    throw new Error('Unable to exec command. Container is not running');
  }

  async isUpdateScriptNeeded(): Promise<boolean> {
    const result = await this.execDockerCommand(['/update_check.sh']);
    if (!result.error) {
      return !!result.updateNeeded;
    }
    throw new Error(result.message);
  }

  async generateAuthToken(): Promise<string> {
    const result = await this.execDockerCommand(['/genkey.sh', '--json']);
    if (!result.error) {
      return result.data;
    }
    throw new Error(result.message);
  }

  async listPackages(): Promise<Package[]> {
    const result = await this.execDockerCommand(
      ['php', '/rnadetector/ws/artisan', 'packages:list'],
      1000
    );
    if (!result.error) {
      return result.packages;
    }
    throw new Error(result.message);
  }

  async execDockerCommandLive(
    Cmd: string[],
    outputCallback: string => void,
    errCallback: ?(string) => void = null,
    exitCallback: ?(number) => void = null,
    timeoutRunning: number = 30000
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
        onExit,
        () => {
          return new Promise(resolve => {
            exec.inspect((e, d) => {
              resolve(d && d.Running);
            });
          });
        },
        timeoutRunning
      );
    }
    throw new Error('Unable to exec command. Container is not running');
  }

  async clearQueue(): Promise<*> {
    const status = await this.checkContainerStatus();
    if (status === 'running') {
      return this.execDockerCommand(
        ['php', '/rnadetector/ws/artisan', 'queue:clear'],
        1000,
        false
      );
    }
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
            status.pushEvent(event);
            if (outputCallback) {
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
  }
  return instance;
};

export const resetInstance = () => {
  instance = null;
};

export default getInstance();
