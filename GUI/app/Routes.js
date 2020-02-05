import React from 'react';
import { Switch, Route } from 'react-router';
import routes from './constants/routes';
import App from './containers/App';
import HomePage from './containers/HomePage';
import CounterPage from './containers/CounterPage';
import SettingsPage from './containers/SettingsPage';
import JobsListPage from './containers/JobsListPage';
import AnnotationsListPage from './containers/AnnotationsListPage';
import ReferencesListPage from './containers/ReferencesListPage';
import CreateReferencePage from './containers/CreateReferencePage';
import CreateAnnotationPage from './containers/CreateAnnotationPage';
import LongRNA from './containers/Analysis/LongRNAPage';
import SmallRNA from './containers/Analysis/SmallRNAPage';
import CircRNA from './containers/Analysis/CircRNAPage';
import SampleGroup from './containers/Analysis/SampleGroupPage';

export default () => (
  <App>
    <Switch>
      <Route path={routes.CREATE_SAMPLE_GROUP} component={SampleGroup} />
      <Route path={routes.ANALYSIS_CIRC_RNA} component={CircRNA} />
      <Route path={routes.ANALYSIS_SMALL_RNA} component={SmallRNA} />
      <Route path={routes.ANALYSIS_LONG_RNA} component={LongRNA} />
      <Route path={routes.SETTINGS} component={SettingsPage} />
      <Route path={routes.CREATE_REFERENCE} component={CreateReferencePage} />
      <Route path={routes.REFERENCES} component={ReferencesListPage} />
      <Route path={routes.JOBS} component={JobsListPage} />
      <Route path={routes.CREATE_ANNOTATION} component={CreateAnnotationPage} />
      <Route path={routes.ANNOTATIONS} component={AnnotationsListPage} />
      <Route path={routes.COUNTER} component={CounterPage} />
      <Route path={routes.HOME} component={HomePage} />
    </Switch>
  </App>
);
