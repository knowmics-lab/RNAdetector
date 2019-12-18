// @flow
import * as JobsActions from '../actions/jobs';
import { Action } from './types';
import * as Api from '../api';
import type { JobsListType } from './types';

const initJobsListState = (): JobsListType => ({
  jobsList: {
    state: {
      current_page: null,
      last_page: null,
      per_page: 15,
      total: null,
      isFetching: false,
      isError: false,
      errorMessage: ''
    },
    pages: {}
  }
});

export default function jobs(
  state: JobsListType,
  action: Action
): JobsListType {
  const oldState = state || initJobsListState();
  switch (action.type) {
    default:
      return oldState;
  }
}
