/* eslint-disable camelcase */
// @flow
import { has } from 'lodash';
import type { Action, Dispatch, GetState } from '../reducers/types';
import type { Job, JobsCollection } from '../types/jobs';
import * as Api from '../api';
import { pushNotificationSimple } from './notifications';
import type { SortingSpec } from '../types/common';

export const JOBS_LIST_LOADING = 'JOBS--LIST--LOADING';
export const JOBS_LIST_LOADED = 'JOBS--LIST--LOADED';
export const JOBS_LIST_CACHED = 'JOBS--LIST--CACHED';
export const JOBS_LIST_ERROR = 'JOBS--LIST--ERROR';
export const JOBS_LIST_SET_PER_PAGE = 'JOBS--LIST--SET_PER_PAGE';
export const JOBS_LIST_SET_SORTING = 'JOBS--LIST--SET_SORTING';
export const JOBS_LIST_RESET_ALL = 'JOBS--LIST--RESET_ALL';
export const JOBS_LIST_RESET_SELECTED = 'JOBS--LIST--RESET_SELECTED';
export const JOBS_LIST_REQUEST_REFRESH = 'JOBS--LIST--REQUEST_REFRESH';
export const JOB_LOADING = 'JOBS--JOB-LOADING';
export const JOB_SUBMITTING = 'JOBS--JOB-SUBMITTING';
export const JOB_SUBMITTED = 'JOBS--JOB-SUBMITTED';
export const JOB_LOADED = 'JOBS--JOB-LOADED';
export const JOB_CACHED = 'JOBS--JOB-CACHED';
export const JOB_ERROR = 'JOBS--JOB-ERROR';
export const JOB_DELETING = 'JOBS--JOB-DELETING';
export const JOB_UPDATE_DELETING = 'JOBS--JOB-UPDATE-DELETING';

let waitForDeleteTimer = null;

export function setPerPage(perPage: number = 15) {
  return async (dispatch: Dispatch, getState: GetState) => {
    const {
      jobs: {
        jobsList: { per_page: oldPerPage }
      }
    } = getState();
    if (oldPerPage !== perPage) {
      dispatch(internalSetPerPage(perPage));
      dispatch(requestPage(1));
    }
  };
}

export function setSorting(sorting: SortingSpec = { created_at: 'desc' }) {
  return async (dispatch: Dispatch) => {
    dispatch(internalSetSorting(sorting));
    dispatch(requestPage(1));
  };
}

export function requestPage(page: number) {
  return async (dispatch: Dispatch, getState: GetState) => {
    try {
      const {
        jobs: {
          jobsList: { refreshAll, refreshPages }
        }
      } = getState();
      if (refreshAll) dispatch(jobsListResetAll());
      if (refreshPages.length > 0)
        dispatch(jobsListResetSelected(refreshPages));
      dispatch(jobsListLoading());
      const {
        jobs: {
          jobsList: {
            pages,
            state: { per_page, sorting }
          }
        }
      } = getState();
      if (!has(pages, page)) {
        const jobs = await Api.Jobs.fetchJobs(per_page, sorting, page);
        dispatch(jobsListLoaded(jobs));
      } else {
        dispatch(jobsListCached(page));
      }
    } catch (e) {
      dispatch(jobsListError());
      dispatch(
        pushNotificationSimple(`An error occurred: ${e.message}!`, 'error')
      );
    }
  };
}

export function refreshPage(page: ?number = null) {
  return (dispatch: Dispatch) => {
    dispatch(jobsListRequestRefresh(!page ? null : [page]));
    if (page) dispatch(requestPage(page));
  };
}

export function requestJob(jobId: number, force: boolean = false) {
  return async (dispatch: Dispatch, getState: GetState) => {
    try {
      dispatch(jobLoading());
      const {
        jobs: {
          jobs: { items }
        }
      } = getState();
      if (!has(items, jobId) || force) {
        const job = await Api.Jobs.fetchJobById(jobId);
        dispatch(jobLoaded(job));
      } else {
        dispatch(jobCached());
      }
    } catch (e) {
      dispatch(jobError());
      dispatch(
        pushNotificationSimple(`An error occurred: ${e.message}!`, 'error')
      );
    }
  };
}

export function submitJob(jobId: number, page: ?number = null) {
  return async (dispatch: Dispatch) => {
    try {
      dispatch(jobSubmitting(jobId));
      const job = await Api.Jobs.submitJob(jobId);
      dispatch(jobLoaded(job));
      if (page) {
        dispatch(refreshPage(page));
      }
      dispatch(pushNotificationSimple(`Job "${job.name}" has been submitted!`));
    } catch (e) {
      dispatch(jobError());
      dispatch(
        pushNotificationSimple(`An error occurred: ${e.message}!`, 'error')
      );
    } finally {
      dispatch(jobSubmitted(jobId));
    }
  };
}

export function waitForDelete() {
  return async (dispatch: Dispatch, getState: GetState) => {
    if (waitForDeleteTimer !== null) return;
    waitForDeleteTimer = setInterval(() => {
      const state = getState();
      const {
        jobs: {
          jobs: { deleting }
        }
      } = state;
      if (deleting.length === 0) {
        clearInterval(waitForDeleteTimer);
      } else {
        // eslint-disable-next-line promise/catch-or-return
        Api.Jobs.processDeletedList(deleting).then(d => {
          // eslint-disable-next-line promise/always-return
          if (d.length === 0) clearInterval(waitForDeleteTimer);
          dispatch(jobUpdateDeleting(d));
          dispatch(refreshPage());
        });
      }
    }, 10000);
  };
}

export function deleteJob(jobId: number) {
  return async (dispatch: Dispatch) => {
    try {
      await Api.Jobs.deleteJob(jobId);
      dispatch(pushNotificationSimple('Job cancellation request sent.'));
      dispatch(jobDeleting(jobId));
      dispatch(waitForDelete());
    } catch (e) {
      dispatch(
        pushNotificationSimple(`An error occurred: ${e.message}!`, 'error')
      );
    }
  };
}

function internalSetPerPage(perPage: number): Action {
  return {
    type: JOBS_LIST_SET_PER_PAGE,
    payload: {
      per_page: perPage
    }
  };
}

function internalSetSorting(payload: SortingSpec): Action {
  return {
    type: JOBS_LIST_SET_SORTING,
    payload
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

export function jobsListError(): Action {
  return {
    type: JOBS_LIST_ERROR,
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

export function jobLoading(): Action {
  return {
    type: JOB_LOADING,
    payload: {}
  };
}

export function jobSubmitting(payload: number): Action {
  return {
    type: JOB_SUBMITTING,
    payload
  };
}

export function jobSubmitted(payload: number): Action {
  return {
    type: JOB_SUBMITTED,
    payload
  };
}

export function jobDeleting(payload: number): Action {
  return {
    type: JOB_DELETING,
    payload
  };
}

export function jobUpdateDeleting(payload: number[]): Action {
  return {
    type: JOB_UPDATE_DELETING,
    payload
  };
}

export function jobLoaded(payload: Job): Action {
  return {
    type: JOB_LOADED,
    payload
  };
}

export function jobCached(): Action {
  return {
    type: JOB_CACHED,
    payload: {}
  };
}

export function jobError(): Action {
  return {
    type: JOB_ERROR,
    payload: {}
  };
}
