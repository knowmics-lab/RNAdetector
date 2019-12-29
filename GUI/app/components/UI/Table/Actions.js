// @flow
import * as React from 'react';
import type { ReadOnlyData, RowAction } from './types';
import Action from './Action';

type Props = {
  actions: RowAction[],
  data: ReadOnlyData,
  keyField: string,
  size: string
};

export default function Actions({ actions, data, keyField, size }: Props) {
  const k = i => `action-${data[keyField]}-${i}`;
  return (
    <>
      {actions.map((a, i) => (
        <Action action={a} data={data} size={size} key={k(i)} />
      ))}
    </>
  );
}
