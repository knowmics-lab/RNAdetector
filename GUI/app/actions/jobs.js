/* eslint-disable camelcase */
// @flow
import { has } from 'lodash';
import type { Action, Dispatch, GetState } from '../reducers/types';
import type { JobsCollection } from '../types/jobs';
import * as Api from '../api';

export const JOBS_LIST_LOADING = 'JOBS--LIST--LOADING';
export const JOBS_LIST_LOADED = 'JOBS--LIST--LOADED';
export const JOBS_LIST_CACHED = 'JOBS--LIST--CACHED';
export const JOBS_LIST_ERROR = 'JOBS--LIST--ERROR';
export const JOBS_LIST_RESET_LOADING = 'JOBS--LIST--RESET_LOADING';
export const JOBS_LIST_SET_PER_PAGE = 'JOBS--LIST--SET_PER_PAGE';
export const JOBS_LIST_RESET_ALL = 'JOBS--LIST--RESET_ALL';
export const JOBS_LIST_RESET_SELECTED = 'JOBS--LIST--RESET_SELECTED';
export const JOBS_LIST_REQUEST_REFRESH = 'JOBS--LIST--REQUEST_REFRESH';

export function setPerPage(perPage: number = 15) {
  return async (dispatch: Dispatch, getState: GetState) => {
    const {
      jobs: {
        jobsList: { per_page: oldPerPage }
      }
    } = getState();
    if (oldPerPage !== perPage) {
      dispatch(internalSetPerPage(perPage));
      dispatch(jobsListRequestRefresh());
      dispatch(requestPage(1));
    }
  };
}

export function requestPage(page: number) {
  return async (dispatch: Dispatch, getState: GetState) => {
    try {
      dispatch(jobsListLoading());
      const {
        jobs: {
          jobsList: { refreshAll, refreshPages }
        }
      } = getState();
      if (refreshAll) dispatch(jobsListResetAll());
      if (refreshPages.length > 0)
        dispatch(jobsListResetSelected(refreshPages));
      const {
        jobs: {
          jobsList: {
            pages,
            state: { per_page }
          }
        }
      } = getState();
      if (!has(pages, page)) {
        const jobs = await Api.Jobs.fetchJobs(per_page, page);
        dispatch(jobsListLoaded(jobs));
      } else {
        dispatch(jobsListCached(page));
      }
    } catch (e) {
      dispatch(jobsListError(e.message));
    }
  };
}

export function internalSetPerPage(perPage: number): Action {
  return {
    type: JOBS_LIST_SET_PER_PAGE,
    payload: {
      per_page: perPage
    }
  };
}

export function jobsListLoading(): Action {
  return {
    type: JOBS_LIST_LOADING,
    payload: {}
  };
}

export function jobsListLoaded(payload: JobsCollection): Action {
  return {
    type: JOBS_LIST_LOADED,
    payload
  };
}

export function jobsListCached(page: number): Action {
  return {
    type: JOBS_LIST_CACHED,
    payload: {
      page
    }
  };
}

export function jobsListError(message: string): Action {
  return {
    type: JOBS_LIST_ERROR,
    payload: {
      message
    }
  };
}

export function jobsListResetLoading(): Action {
  return {
    type: JOBS_LIST_RESET_LOADING,
    payload: {}
  };
}

export function jobsListResetAll(): Action {
  return {
    type: JOBS_LIST_RESET_ALL,
    payload: {}
  };
}

export function jobsListResetSelected(pages: number[]): Action {
  return {
    type: JOBS_LIST_RESET_SELECTED,
    payload: { pages }
  };
}

export function jobsListRequestRefresh(pages: ?(number[]) = null): Action {
  return {
    type: JOBS_LIST_REQUEST_REFRESH,
    payload: {
      all: pages !== null,
      pages
    }
  };
}
