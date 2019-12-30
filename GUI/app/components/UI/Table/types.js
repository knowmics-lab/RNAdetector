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

export type TableState = {
  currentPage: ?number,
  rowsPerPage: ?number,
  totalRows: ?number,
  isLoading: boolean
};

export type RowActionFunction = ReadOnlyData => ChildrenArray<ReactElement<*>>;

export type RowActionObject = {
  shown: boolean | (ReadOnlyData => boolean),
  icon: string | (() => ReactElement<*>),
  size?: string,
  disabled?: boolean | (ReadOnlyData => boolean),
  color?: string,
  onClick?: (MouseEvent, ReadOnlyData) => void,
  tooltip?: string
};

export type RowActionType = RowActionObject | RowActionFunction;

export type ToolbarActionFunction = TableState => ChildrenArray<
  ReactElement<*>
>;

export type ToolbarActionCustom = {
  custom: true,
  action: ToolbarActionFunction,
  align: 'left' | 'center' | 'right'
};

export type ToolbarActionButton = {
  custom?: false,
  align: 'left' | 'center' | 'right',
  shown: boolean | (TableState => boolean),
  icon: string | (() => ReactElement<*>),
  disabled?: boolean | (TableState => boolean),
  color?: string,
  onClick?: (MouseEvent, TableState) => void,
  tooltip?: string
};

export type ToolbarActionType = ToolbarActionCustom | ToolbarActionButton;
