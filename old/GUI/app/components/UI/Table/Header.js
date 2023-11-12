/* eslint-disable no-nested-ternary */
// @flow
import React from 'react';
import { has } from 'lodash';
import { makeStyles } from '@material-ui/core/styles';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import TableSortLabel from '@material-ui/core/TableSortLabel';
import TableCell from '@material-ui/core/TableCell';
import Checkbox from '@material-ui/core/Checkbox';
import type { TableColumn } from './types';
import type { SortingSpec } from '../../../types/common';
import * as Api from '../../../api';

const useStyles = makeStyles({
  stickyStyle: {
    backgroundColor: 'white'
  }
});

type Props = {
  columns: TableColumn[],
  sorting: SortingSpec,
  sortable: boolean,
  changeSorting: SortingSpec => void,
  hasCheckbox?: boolean,
  single?: boolean,
  selectedAny?: boolean,
  selectedAll?: boolean,
  handleSelect?: () => void
};

export default function Header({
  columns,
  sorting,
  sortable,
  changeSorting,
  hasCheckbox,
  single,
  selectedAny,
  selectedAll,
  handleSelect
}: Props) {
  const classes = useStyles();
  const sf = column => column.sortingField || column.dataField;
  const makeChangeHandler = column => event => {
    const hasSorting = has(sorting, column);
    const oldDirection = hasSorting ? sorting[column] : null;
    const newDirection =
      oldDirection === null ? 'desc' : oldDirection === 'desc' ? 'asc' : null;
    const newSorting =
      newDirection === null
        ? Api.Utils.filterByKey(sorting, k => k !== column)
        : { ...sorting, [column]: newDirection };
    changeSorting(newSorting);
    event.preventDefault();
  };
  return (
    <TableHead>
      <TableRow>
        {hasCheckbox && (
          <TableCell padding="checkbox">
            {!single && (
              <Checkbox
                indeterminate={selectedAny}
                checked={selectedAll}
                onChange={handleSelect}
              />
            )}
          </TableCell>
        )}
        {columns.map(column =>
          column !== 'actions' ? (
            <TableCell
              key={column.dataField}
              align={column.align}
              style={{ minWidth: column.minWidth }}
              className={classes.stickyStyle}
            >
              {sortable && !column.disableSorting ? (
                <TableSortLabel
                  active={has(sorting, sf(column))}
                  direction={
                    has(sorting, sf(column)) ? sorting[sf(column)] : 'desc'
                  }
                  onClick={makeChangeHandler(sf(column))}
                >
                  {column.label}
                </TableSortLabel>
              ) : (
                column.label
              )}
            </TableCell>
          ) : (
            <TableCell
              key="actions"
              align="center"
              className={classes.stickyStyle}
            >
              Actions
            </TableCell>
          )
        )}
      </TableRow>
    </TableHead>
  );
}

Header.defaultProps = {
  hasCheckbox: false,
  selectedAny: false,
  selectedAll: false,
  handleSelect: undefined
};
