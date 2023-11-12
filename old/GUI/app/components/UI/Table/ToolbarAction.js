// @flow
import * as React from 'react';
import Icon from '@material-ui/core/Icon';
import IconButton from '../IconButton';
import type { TableState, ToolbarActionType } from './types';

export type Props = {
  action: ToolbarActionType,
  state: TableState
};

function isF(x: mixed): boolean %checks {
  return typeof x === 'function';
}

export default function ToolbarAction({ action, state }: Props) {
  if (action.custom) {
    if (action.action && isF(action.action)) {
      return action.action(state);
    }
    return null;
  }
  const shown = isF(action.shown) ? action.shown(state) : action.shown;
  if (!shown) return null;
  const actionDisabled = action.disabled || false;
  const disabled = isF(actionDisabled) ? actionDisabled(state) : actionDisabled;
  const icon = isF(action.icon) ? (
    action.icon()
  ) : (
    <Icon className={action.icon} fontSize="inherit" />
  );
  const color = action.color || 'inherit';
  const onClick = event =>
    action.onClick ? action.onClick(event, state) : undefined;
  return (
    <IconButton
      color={color}
      disabled={disabled}
      onClick={onClick}
      title={action.tooltip}
    >
      {icon}
    </IconButton>
  );
}
