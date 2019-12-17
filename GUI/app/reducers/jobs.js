// @flow
import * as JobsActions from '../actions/jobs';
import { Action } from './types';
import * as Api from '../api';
import type { jobsListType } from './types';

const initJobsListState = (): jobsListType => ({
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
  state: ?jobsListType = initJobsListState(),
  action: Action
) {
  switch (action.type) {
    default:
      return state;
  }
}
