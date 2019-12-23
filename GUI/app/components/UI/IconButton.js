// @flow
import React, { useMemo, forwardRef } from 'react';
import { Link as RouterLink } from 'react-router-dom';
import { Tooltip, IconButton as IB } from '@material-ui/core';

export type IconButtonType = {
  color?: 'default' | 'inherit' | 'primary' | 'secondary',
  title: string,
  children: *,
  disabled?: boolean,
  href?: ?string,
  size?: 'small' | 'medium',
  onClick?: ?() => void
};

export default function IconButton({
  color,
  title,
  children,
  disabled,
  href,
  size,
  onClick
}: IconButtonType) {
  if (!href && onClick) {
    return (
      <Tooltip title={title}>
        <IB color={color} disabled={disabled} onClick={onClick} size={size}>
          {children}
        </IB>
      </Tooltip>
    );
  }
  if (!onClick && href) {
    const renderLink = useMemo(
      () =>
        forwardRef((itemProps, ref) => (
          <RouterLink
            to={href}
            // eslint-disable-next-line react/jsx-props-no-spreading
            {...itemProps}
            innerRef={ref}
          />
        )),
      [href]
    );
    return (
      <Tooltip title={title}>
        <IB
          color={color}
          disabled={disabled}
          component={renderLink}
          size={size}
        >
          {children}
        </IB>
      </Tooltip>
    );
  }
  throw new Error('You must specify onClick or href.');
}

IconButton.defaultProps = {
  color: 'default',
  size: 'medium',
  variant: 'text',
  disabled: false,
  href: null,
  onClick: null
};
