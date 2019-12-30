// @flow
import * as React from 'react';
import type { ReadOnlyData, RowActionType } from './types';
import RowAction from './RowAction';

type Props = {
  actions: RowActionType[],
  data: ReadOnlyData,
  keyField: string,
  size: string
};

export default function RowActions({ actions, data, keyField, size }: Props) {
  const k = i => `action-${data[keyField]}-${i}`;
  return (
    <>
      {actions.map((a, i) => (
        <RowAction action={a} data={data} size={size} key={k(i)} />
      ))}
    </>
  );
}
