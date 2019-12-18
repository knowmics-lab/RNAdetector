// @flow
import type { Action, Dispatch } from '../reducers/types';
import * as Api from '../api';
import type { ConfigObjectType } from '../types/settings';

export const SETTINGS_SAVING = 'SETTINGS--SAVING';
export const SETTINGS_SAVED = 'SETTINGS--SAVED';
export const SETTINGS_ERROR = 'SETTINGS--ERROR';
export const SETTINGS_RESET_SAVED = 'SETTINGS--RESET-SAVED';

export function saveSettings(newSettings: ConfigObjectType): Action {
  return async (dispatch: Dispatch) => {
    try {
      dispatch(settingsSaving());
      dispatch(settingsSaved(await Api.Settings.saveConfig(newSettings)));
    } catch (e) {
      dispatch(settingsError(e));
    }
  };
}

export function settingsSaving(): Action {
  return {
    type: SETTINGS_SAVING,
    payload: {}
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
