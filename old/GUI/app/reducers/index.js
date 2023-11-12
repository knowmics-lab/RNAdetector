// @flow
import { combineReducers } from 'redux';
import { connectRouter } from 'connected-react-router';
import type { HashHistory } from 'history';
import counter from './counter';
import annotations from './annotations';
import jobs from './jobs';
import references from './references';
import settings from './settings';
import notifications from './notifications';

export default function createRootReducer(history: HashHistory) {
  return combineReducers<{}, *>({
    router: connectRouter(history),
    counter,
    settings,
    annotations,
    jobs,
    references,
    notifications
  });
}
