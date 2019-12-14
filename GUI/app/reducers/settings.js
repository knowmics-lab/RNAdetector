// @flow
import {
  SETTINGS_SAVED,
  SETTINGS_ERROR,
  SETTINGS_RESET_SAVED
} from '../actions/settings';
import { settingsStateType, Action } from './types';
import * as Api from '../api';

const initConfigState = (): settingsStateType => ({
  state: {
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
    case SETTINGS_ERROR:
      return {
        ...state,
        state: {
          saved: false,
          error: true,
          message: action.payload.message
        }
      };
    case SETTINGS_SAVED:
      return {
        state: {
          ...state.state,
          saved: true
        },
        ...action.payload.settings
      };
    case SETTINGS_RESET_SAVED:
      return {
        ...state,
        state: {
          saved: false,
          error: false,
          message: ''
        }
      };
    default:
      return state;
  }
}
