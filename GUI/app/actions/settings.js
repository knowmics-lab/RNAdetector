// @flow
import type { Action, Dispatch } from '../reducers/types';
import type { ConfigObjectType } from '../api';
import * as Api from '../api';

export const SETTINGS_SAVED = 'SETTINGS--SAVED';
export const SETTINGS_ERROR = 'SETTINGS--ERROR';
export const SETTINGS_RESET_SAVED = 'SETTINGS--RESET-SAVED';

export function saveSettings(newSettings: ConfigObjectType): Action {
  return async (dispatch: Dispatch) => {
    try {
      dispatch(settingsSaved(await Api.Settings.saveConfig(newSettings)));
    } catch (e) {
      dispatch(settingsError(e));
    }
  };
}

export function settingsSaved(newSettings: ConfigObjectType): Action {
  return {
    type: SETTINGS_SAVED,
    payload: {
      settings: newSettings
    }
  };
}

export function settingsError(message: string): Action {
  return {
    type: SETTINGS_ERROR,
    payload: {
      message
    }
  };
}

export function resetSaved(): Action {
  return {
    type: SETTINGS_RESET_SAVED,
    payload: {}
  };
}
