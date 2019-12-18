// @flow
import {
  SETTINGS_SAVING,
  SETTINGS_SAVED,
  SETTINGS_ERROR,
  SETTINGS_RESET_SAVED
} from '../actions/settings';
import { Action } from './types';
import * as Api from '../api';
import type { SettingsStateType } from '../types/settings';

const initConfigState = (): SettingsStateType => ({
  state: {
    saving: false,
    saved: false,
    error: false,
    message: ''
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
          saving: true,
          saved: false,
          error: false,
          message: ''
        }
      };
    case SETTINGS_ERROR:
      return {
        ...oldState,
        state: {
          saving: false,
          saved: false,
          error: true,
          message: action.payload.message
        }
      };
    case SETTINGS_SAVED:
      return {
        state: {
          saving: false,
          saved: true,
          error: false,
          message: ''
        },
        ...action.payload.settings
      };
    case SETTINGS_RESET_SAVED:
      return {
        ...oldState,
        state: {
          saving: false,
          saved: false,
          error: false,
          message: ''
        }
      };
    default:
      return oldState;
  }
}
