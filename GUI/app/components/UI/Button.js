// @flow
import React, { useMemo, forwardRef } from 'react';
import { Link as RouterLink } from 'react-router-dom';
import { Button as OB } from '@material-ui/core';

export type ButtonType = {
  color?: 'default' | 'inherit' | 'primary' | 'secondary',
  children: *,
  disabled?: boolean,
  href?: ?string,
  size?: 'small' | 'medium' | 'large',
  variant?: 'text' | 'outlined' | 'contained',
  onClick?: ?() => void
};

export default function Button({
  color: co,
  children: c,
  disabled: d,
  href,
  size: s,
  variant: v,
  onClick
}: ButtonType) {
  if (!href && onClick) {
    return (
      <OB variant={v} color={co} disabled={d} onClick={onClick} size={s}>
        {c}
      </OB>
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
      <OB variant={v} color={co} disabled={d} size={s} component={renderLink}>
        {c}
      </OB>
    );
  }
  throw new Error('You must specify onClick or href.');
}

Button.defaultProps = {
  color: 'default',
  size: 'medium',
  variant: 'text',
  disabled: false,
  href: null,
  onClick: null
};
