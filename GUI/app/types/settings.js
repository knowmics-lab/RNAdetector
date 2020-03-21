// @flow
import type { ModifiableStateType } from './common';

export type ConfigObjectType = {
  +configured?: boolean,
  +local: boolean,
  +apiProtocol: 'http' | 'https',
  +apiHostname: string,
  +apiPort: number,
  +apiPath: string,
  +publicPath: string,
  +dataPath: string,
  +containerName: string,
  +apiKey: string,
  +dockerExecutablePath: string,
  +socketPath?: ?string
};

export type SettingsStateType = {
  +settings: {
    +state: ModifiableStateType
  } & ConfigObjectType
};
