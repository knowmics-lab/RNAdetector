/* eslint-disable no-case-declarations */
// @flow
import * as JobsActions from '../actions/jobs';
import { Action } from './types';
import * as Api from '../api';
import type { JobsListType, JobsStateType } from './types';

const initJobsListState = (perPage: number = 15): JobsStateType => ({
  jobsList: {
    refreshAll: false,
    refreshPages: [],
    state: {
      current_page: null,
      last_page: null,
      per_page: perPage,
      total: null,
      isFetching: false,
      isError: false,
      errorMessage: ''
    },
    pages: {}
  }
});

function resetSelected(state: JobsListType, pages: number[]): JobsListType {
  return {
    ...state,
    refreshPages: [],
    pages: Api.Utils.filterByKey(state.pages, key => pages.includes(key))
  };
}

function jobsList(state: JobsListType, action: Action): JobsListType {
  switch (action.type) {
    case JobsActions.JOBS_LIST_RESET_LOADING:
      return {
        ...state,
        state: {
          ...state.state,
          isFetching: false,
          isError: false,
          errorMessage: ''
        }
      };
    case JobsActions.JOBS_LIST_RESET_ALL:
      return initJobsListState(state.state.per_page);
    case JobsActions.JOBS_LIST_RESET_SELECTED:
      return resetSelected(state, action.payload.pages);
    case JobsActions.JOBS_LIST_ERROR:
      const {
        payload: { message }
      } = action;
      return {
        ...state,
        state: {
          ...state.state,
          isFetching: false,
          isError: true,
          errorMessage: message
        }
      };
    case JobsActions.JOBS_LIST_SET_PER_PAGE:
      return {
        ...state,
        refreshAll: true,
        refreshPages: [],
        state: {
          ...state.state,
          per_page: action.payload.per_page
        }
      };
    default:
      return state;
  }
}

export default function jobs(
  state: JobsStateType,
  action: Action
): JobsStateType {
  const oldState = state || initJobsListState();
  return {
    jobsList: jobsList(oldState.jobsList, action)
  };
}
