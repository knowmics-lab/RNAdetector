/* eslint-disable react/no-array-index-key */
// @flow
import * as React from 'react';
import MaterialToolbar from '@material-ui/core/Toolbar';
import { makeStyles } from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import type { TableState, ToolbarActionType } from './types';
import ToolbarAction from './ToolbarAction';

const useStyles = makeStyles(theme => ({
  root: {
    paddingRight: theme.spacing(1),
    flexFlow: 'row nowrap',
    justifyContent: 'space-between',
    alignItems: 'center'
  },
  title: {
    flex: '0 0 auto',
    color: theme.palette.text.primary,
    order: 0
  },
  left: {
    flex: '1 1 auto',
    color: theme.palette.text.secondary,
    order: 1
  },
  right: {
    flex: '1 1 auto',
    color: theme.palette.text.secondary,
    textAlign: 'right',
    order: 3
  },
  center: {
    flex: '1 0 auto',
    color: theme.palette.text.secondary,
    textAlign: 'center',
    order: 2
  },
  actions: {
    color: theme.palette.text.secondary
  }
}));

type Props = {
  title?: ?(string | React.ChildrenArray<React.Element<*>>),
  actions: ToolbarActionType[],
  state: TableState
};

Toolbar.defaultProps = {
  title: null
};

export default function Toolbar({ title, actions, state }: Props) {
  const classes = useStyles();
  if (!actions || actions.length === 0) return null;
  const renderTitle = () => {
    if (!title) return null;
    return (
      <div className={classes.title}>
        {typeof title === 'string' ? (
          <Typography variant="h6">{title}</Typography>
        ) : (
          title
        )}
      </div>
    );
  };
  const renderActions = d =>
    actions
      .filter(a => a.align === d)
      .map((a, i) => (
        <ToolbarAction action={a} state={state} key={`toolbar-${d}-${i}`} />
      ));
  return (
    <MaterialToolbar className={classes.root}>
      {renderTitle()}
      <div className={classes.left}>{renderActions('left')}</div>
      <div className={classes.center}>{renderActions('center')}</div>
      <div className={classes.right}>{renderActions('right')}</div>
    </MaterialToolbar>
  );
}
