// @flow
import React from 'react';
import Icon from '@material-ui/core/Icon';
import IconButton from '../IconButton';
import type { ReadOnlyData, RowActionType } from './types';

export type Props = {
  action: RowActionType,
  data: ReadOnlyData,
  size: string
};

function isF(x: mixed): boolean %checks {
  return typeof x === 'function';
}

export default function RowAction({ action, data, size }: Props) {
  if (typeof action === 'function') {
    return <>{action(data, size)}</>;
  }
  const shown = isF(action.shown) ? action.shown(data) : action.shown;
  if (!shown) return null;
  const actionDisabled = action.disabled || false;
  const disabled = isF(actionDisabled) ? actionDisabled(data) : actionDisabled;
  const icon = isF(action.icon) ? (
    action.icon()
  ) : (
    <Icon className={action.icon} fontSize="inherit" />
  );
  const color = action.color || 'inherit';
  const onClick = event =>
    action.onClick ? action.onClick(event, data) : undefined;
  return (
    <IconButton
      size={action.size || size}
      color={color}
      disabled={disabled}
      onClick={onClick}
      title={action.tooltip}
    >
      {icon}
    </IconButton>
  );
}
