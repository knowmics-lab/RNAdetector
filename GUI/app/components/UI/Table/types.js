// @flow
import type { ChildrenArray, Element as ReactElement } from 'react';

export type ReadOnlyData = { +[string]: * };

export type TableColumn =
  | 'actions'
  | {
      id: string,
      label: string,
      minWidth?: number,
      align?: 'left' | 'right' | 'center' | 'justify',
      format?: ReadOnlyData => mixed
    };

export type RowActionFunction = ReadOnlyData => ChildrenArray<ReactElement<*>>;

export type RowActionObject = {
  shown: boolean | (ReadOnlyData => boolean),
  icon: string | (() => ReactElement<*>),
  size?: string,
  disabled?: boolean | (ReadOnlyData => boolean),
  color?: string,
  onClick?: (MouseEvent, ReadOnlyData) => void,
  href?: string,
  tooltip?: string
};

export type RowAction = RowActionObject | RowActionFunction;
