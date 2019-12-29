// @flow
import React, { useMemo, forwardRef } from 'react';
import type { ChildrenArray, Element } from 'react';
import { Link as RouterLink } from 'react-router-dom';
import { Tooltip, IconButton as IB } from '@material-ui/core';

export type IconButtonType = {
  color?: string,
  title?: string,
  children: ChildrenArray<Element<*>>,
  disabled?: boolean,
  href?: ?string,
  size?: string,
  onClick?: ?(MouseEvent) => void
};

const makeTooltip = (component, tooltip, disabled) => {
  if (!tooltip) return component;
  return disabled ? (
    <Tooltip title={tooltip}>
      <span>{component}</span>
    </Tooltip>
  ) : (
    <Tooltip title={tooltip}>{component}</Tooltip>
  );
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
    return makeTooltip(
      <IB color={color} disabled={disabled} component={renderLink} size={size}>
        {children}
      </IB>,
      title,
      disabled
    );
  }
  const onClickFn = onClick || (() => undefined);
  return makeTooltip(
    <IB color={color} disabled={disabled} onClick={onClickFn} size={size}>
      {children}
    </IB>,
    title,
    disabled
  );
}

IconButton.defaultProps = {
  color: 'default',
  size: 'medium',
  disabled: false,
  title: '',
  href: null,
  onClick: null
};
