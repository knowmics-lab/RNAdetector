/* eslint-disable no-case-declarations */
// @flow
import * as JobsActions from '../actions/jobs';
import { Action } from './types';
import * as Api from '../api';
import type {
  JobsListType,
  JobsStateType,
  JobsCollection,
  LoadedJobs
} from '../types/jobs';
import type { SortingSpec } from '../types/common';

const initJobsListState = (
  perPage: number = 15,
  sorting: SortingSpec = { created_at: 'desc' }
): JobsStateType => ({
  jobsList: {
    refreshAll: false,
    refreshPages: [],
    state: {
      current_page: null,
      last_page: null,
      per_page: perPage,
      total: null,
      sorting,
      fetching: false
    },
    pages: {}
  },
  jobs: {
    fetching: false,
    submitting: [],
    deleting: [],
    items: {}
  }
});

function resetJobsListSelected(
  state: JobsListType,
  pages: number[]
): JobsListType {
  return {
    ...state,
    refreshPages: [],
    pages: Api.Utils.filterByKey(state.pages, key => !pages.includes(key))
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
      sorting: payload.meta.sorting,
      fetching: false
    },
    pages: {
      ...state.pages,
      [payload.meta.current_page]: Object.values(payload.data)
    }
  };
}

function jobsList(state: JobsListType, action: Action): JobsListType {
  switch (action.type) {
    case JobsActions.JOBS_LIST_RESET_ALL:
      const newState = initJobsListState(
        state.state.per_page,
        state.state.sorting
      );
      return newState.jobsList;
    case JobsActions.JOBS_LIST_RESET_SELECTED:
      return resetJobsListSelected(state, action.payload.pages);
    case JobsActions.JOBS_LIST_ERROR:
      return changeJobsListState(state, {
        current_page: state.state.current_page || 0,
        last_page: state.state.last_page || 0,
        total: state.state.total || 0,
        fetching: false
      });
    case JobsActions.JOBS_LIST_SET_PER_PAGE:
      return {
        ...changeJobsListState(state, { per_page: action.payload.per_page }),
        refreshAll: true,
        refreshPages: []
      };
    case JobsActions.JOBS_LIST_SET_SORTING:
      return {
        ...changeJobsListState(state, { sorting: action.payload }),
        refreshAll: true,
        refreshPages: []
      };
    case JobsActions.JOBS_LIST_LOADING:
      return changeJobsListState(state, {
        fetching: true
      });
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
        fetching: false
      });
    default:
      return state;
  }
}

function changeSubmitting(
  state: LoadedJobs,
  jobId: number,
  isSubmitting: boolean
) {
  let { submitting } = state;
  if (isSubmitting && !submitting.includes(jobId))
    submitting = [...submitting, jobId];
  if (!isSubmitting && submitting.includes(jobId))
    submitting = submitting.filter(k => k !== jobId);
  return {
    ...state,
    submitting
  };
}

function changeFetching(state, fetching) {
  return {
    ...state,
    fetching
  };
}

function loadedJobs(state: LoadedJobs, action: Action): LoadedJobs {
  switch (action.type) {
    case JobsActions.JOB_ERROR:
    case JobsActions.JOB_LOADING:
      return changeFetching(state, true);
    case JobsActions.JOB_SUBMITTING:
      return changeSubmitting(state, action.payload, true);
    case JobsActions.JOB_SUBMITTED:
      return changeSubmitting(state, action.payload, false);
    case JobsActions.JOB_LOADED:
      return {
        fetching: false,
        submitting: state.submitting,
        deleting: state.deleting,
        items: {
          ...state.items,
          [action.payload.id]: action.payload
        }
      };
    case JobsActions.JOB_CACHED:
      return changeFetching(state, false);
    case JobsActions.JOB_DELETING:
      if (state.deleting.includes(action.payload)) return state;
      return {
        ...state,
        deleting: [...state.deleting, action.payload]
      };
    case JobsActions.JOB_UPDATE_DELETING:
      return {
        ...state,
        deleting: action.payload
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
  const newState = {
    jobsList: jobsList(oldState.jobsList, action),
    jobs: loadedJobs(oldState.jobs, action)
  };
  if (
    newState.jobsList === oldState.jobsList &&
    newState.jobs === oldState.jobs
  ) {
    return oldState;
  }
  return newState;
}
