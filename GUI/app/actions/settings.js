// @flow
import type { Action, settingsStateType } from '../reducers/types';

export const SAVE_SETTINGS = 'SAVE_SETTINGS';

export function saveSettings(newSettings: settingsStateType): Action {
  return {
    type: SAVE_SETTINGS,
    data: { settings: newSettings }
  };
}
