// @flow
import {
  SETTINGS_SAVING,
  SETTINGS_SAVED,
  SETTINGS_ERROR
} from '../actions/settings';
import { Action } from './types';
import * as Api from '../api';
import type { SettingsStateType } from '../types/settings';

const initConfigState = (): SettingsStateType => ({
  state: {
    saving: false
  },
  ...Api.Settings.getConfig()
});

export default function settings(
  state: ?SettingsStateType,
  action: Action
): SettingsStateType {
  const oldState = state || initConfigState();
  switch (action.type) {
    case SETTINGS_SAVING:
      return {
        ...oldState,
        state: {
          saving: true
        }
      };
    case SETTINGS_ERROR:
      return {
        ...oldState,
        state: {
          saving: false
        }
      };
    case SETTINGS_SAVED:
      return {
        state: {
          saving: false
        },
        ...action.payload.settings
      };
    default:
      return oldState;
  }
}
