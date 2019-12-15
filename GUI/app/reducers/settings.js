// @flow
import {
  SETTINGS_SAVING,
  SETTINGS_SAVED,
  SETTINGS_ERROR,
  SETTINGS_RESET_SAVED
} from '../actions/settings';
import { settingsStateType, Action } from './types';
import * as Api from '../api';

const initConfigState = (): settingsStateType => ({
  state: {
    saving: false,
    saved: false,
    error: false,
    message: ''
  },
  ...Api.Settings.getConfig()
});

export default function settings(
  state: ?settingsStateType = initConfigState(),
  action: Action
) {
  switch (action.type) {
    case SETTINGS_SAVING:
      return {
        ...state,
        state: {
          saving: true,
          saved: false,
          error: false,
          message: ''
        }
      };
    case SETTINGS_ERROR:
      return {
        ...state,
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
        ...state,
        state: {
          saving: false,
          saved: false,
          error: false,
          message: ''
        }
      };
    default:
      return state;
  }
}
