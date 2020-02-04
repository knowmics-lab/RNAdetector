// @flow
import React from 'react';
import TableBody from '@material-ui/core/TableBody';
import TableRow from '@material-ui/core/TableRow';
import TableCell from '@material-ui/core/TableCell';
import Checkbox from '@material-ui/core/Checkbox';
import type { ReadOnlyData, RowActionType, TableColumn } from './types';
import RowActions from './RowActions';

function Cell(
  column: TableColumn,
  row: ReadOnlyData,
  keyBase: string,
  actions: RowActionType[],
  keyField: string,
  size: string
) {
  if (column !== 'actions') {
    const value = row[column.dataField];
    return (
      <TableCell key={`${keyBase}-${column.dataField}`} align={column.align}>
        {column.format ? column.format(value, row) : value}
      </TableCell>
    );
  }
  return (
    <TableCell key={`${keyBase}-actions`} align="center">
      <RowActions
        actions={actions}
        data={row}
        keyField={keyField}
        size={size}
      />
    </TableCell>
  );
}

type Props = {
  data: ReadOnlyData[],
  keyField: string,
  columns: TableColumn[],
  actions: RowActionType[],
  size: string,
  hasCheckbox?: boolean,
  selectedItems?: boolean[],
  handleSelect?: mixed => void
};

export default function Body({
  data,
  keyField,
  columns,
  actions,
  size,
  hasCheckbox,
  selectedItems,
  handleSelect
}: Props) {
  const isSelected = id =>
    hasCheckbox && selectedItems ? selectedItems.includes(id) : false;
  return (
    <TableBody>
      {Array.isArray(data) &&
        data.length > 0 &&
        data.map(row => {
          const id = row[keyField];
          return (
            <TableRow
              hover
              role="checkbox"
              tabIndex={-1}
              key={`row-${id}`}
              selected={isSelected(id)}
              onClick={() => (handleSelect ? handleSelect(id) : undefined)}
            >
              {hasCheckbox && (
                <TableCell padding="checkbox">
                  <Checkbox checked={isSelected(id)} />
                </TableCell>
              )}
              {columns.map(column =>
                Cell(
                  column,
                  row,
                  `row-${row[keyField]}`,
                  actions,
                  keyField,
                  size
                )
              )}
            </TableRow>
          );
        })}
    </TableBody>
  );
}

Body.defaultProps = {
  hasCheckbox: false,
  selectedItems: [],
  handleSelect: undefined
};
