// @flow
import * as React from 'react';
import { makeStyles, Theme } from '@material-ui/core/styles';
import * as Api from '../../api';
import NotificationsList from '../Layout/Notifications';
import SetupWizard from './SetupWizard';

const useStyles = makeStyles((theme: Theme) => ({
  root: {
    display: 'flex'
  },
  content: {
    flexGrow: 1,
    padding: theme.spacing(3)
  }
}));

type Props = {
  children: React.Node | React.Node[]
};

export default function SetupWizardContainer({ children }: Props) {
  const classes = useStyles();
  return (
    <>
      {Api.Settings.isConfigured() ? (
        children
      ) : (
        <div className={classes.root}>
          <main className={classes.content}>
            <SetupWizard />
            <NotificationsList />
          </main>
        </div>
      )}
    </>
  );
}
