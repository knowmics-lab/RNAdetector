// @flow
import Store from 'electron-store';
import { SAVE_SETTINGS } from '../actions/settings';
import { settingsStateType, Action } from './types';
import configSchema from '../constants/config-schema.json';

const configStore = new Store({ configSchema });

const initConfigState = (): settingsStateType => ({
  local: configStore.get('local', true),
  webserviceUrl: configStore.get('webserviceUrl', 'https://localhost:9898/'),
  jobsPath: configStore.get('jobsPath', '')
});

export default function settings(
  state: ?settingsStateType = initConfigState(),
  action: Action
) {
  switch (action.type) {
    case SAVE_SETTINGS:
      Object.keys(action.data.settings).forEach(prop => {
        configStore.set(prop, action.data.settings[prop]);
      });
      return action.data.settings;
    default:
      return state;
  }
}
