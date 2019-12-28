// @flow
import type { Action, Dispatch } from '../reducers/types';
import * as Api from '../api';
import { pushNotificationSimple } from './notifications';
import type { ConfigObjectType } from '../types/settings';

export const SETTINGS_SAVING = 'SETTINGS--SAVING';
export const SETTINGS_SAVED = 'SETTINGS--SAVED';
export const SETTINGS_ERROR = 'SETTINGS--ERROR';

export function saveSettings(newSettings: ConfigObjectType): Action {
  return async (dispatch: Dispatch) => {
    try {
      dispatch(settingsSaving());
      dispatch(settingsSaved(await Api.Settings.saveConfig(newSettings)));
      dispatch(pushNotificationSimple(`Settings saved!`));
    } catch (e) {
      dispatch(settingsError());
      dispatch(
        pushNotificationSimple(`An error occurred: ${e.message}!`, 'error')
      );
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

export function settingsError(): Action {
  return {
    type: SETTINGS_ERROR,
    payload: {}
  };
}
