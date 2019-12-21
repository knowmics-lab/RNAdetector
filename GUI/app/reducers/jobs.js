/* eslint-disable no-case-declarations */
// @flow
import * as JobsActions from '../actions/jobs';
import { Action } from './types';
import * as Api from '../api';
import type {
  JobsListType,
  JobsStateType,
  JobsCollection
} from '../types/jobs';

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

function resetJobsListSelected(
  state: JobsListType,
  pages: number[]
): JobsListType {
  return {
    ...state,
    refreshPages: [],
    pages: Api.Utils.filterByKey(state.pages, key => pages.includes(key))
  };
}

function changeJobsListState(state, newState) {
  return {
    ...state,
    state: {
      ...state.state,
      ...newState
    }
  };
}

function addLoadedPayload(
  state: JobsListType,
  payload: JobsCollection
): JobsListType {
  return {
    ...state,
    state: {
      ...state.state,
      current_page: payload.meta.current_page,
      last_page: payload.meta.last_page,
      total: payload.meta.total,
      isFetching: false,
      isError: false,
      errorMessage: ''
    },
    pages: {
      ...state.pages,
      [payload.meta.current_page]: Object.values(payload.data)
    }
  };
}

function jobsList(state: JobsListType, action: Action): JobsListType {
  switch (action.type) {
    case JobsActions.JOBS_LIST_RESET_LOADING:
      return changeJobsListState(state, {
        isFetching: false,
        isError: false,
        errorMessage: ''
      });
    case JobsActions.JOBS_LIST_RESET_ALL:
      const newState = initJobsListState(state.state.per_page);
      return newState.jobsList;
    case JobsActions.JOBS_LIST_RESET_SELECTED:
      return resetJobsListSelected(state, action.payload.pages);
    case JobsActions.JOBS_LIST_ERROR:
      return changeJobsListState(state, {
        isFetching: false,
        isError: true,
        errorMessage: action.payload.message
      });
    case JobsActions.JOBS_LIST_SET_PER_PAGE:
      return {
        ...changeJobsListState(state, { per_page: action.payload.per_page }),
        refreshAll: true,
        refreshPages: []
      };
    case JobsActions.JOBS_LIST_LOADING:
      return changeJobsListState(state, { isFetching: true });
    case JobsActions.JOBS_LIST_REQUEST_REFRESH:
      return {
        ...state,
        refreshAll: state.refreshAll || action.payload.all,
        refreshPages: [...state.refreshPages, ...(action.payload.pages || [])]
      };
    case JobsActions.JOBS_LIST_LOADED:
      return addLoadedPayload(state, action.payload);
    case JobsActions.JOBS_LIST_CACHED:
      return changeJobsListState(state, {
        current_page: action.payload.page,
        isFetching: false,
        isError: false,
        errorMessage: ''
      });
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
