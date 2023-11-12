/* eslint-disable camelcase */
// @flow
import { has } from 'lodash';
import type { Action, Dispatch, GetState } from '../reducers/types';
import type { Reference, ReferencesCollection } from '../types/references';
import * as Api from '../api';
import { pushNotificationSimple } from './notifications';
import type { SortingSpec } from '../types/common';

export const REFERENCES_LIST_LOADING = 'REFERENCES--LIST--LOADING';
export const REFERENCES_LIST_LOADED = 'REFERENCES--LIST--LOADED';
export const REFERENCES_LIST_CACHED = 'REFERENCES--LIST--CACHED';
export const REFERENCES_LIST_ERROR = 'REFERENCES--LIST--ERROR';
export const REFERENCES_LIST_SET_PER_PAGE = 'REFERENCES--LIST--SET_PER_PAGE';
export const REFERENCES_LIST_SET_SORTING = 'REFERENCES--LIST--SET_SORTING';
export const REFERENCES_LIST_RESET_ALL = 'REFERENCES--LIST--RESET_ALL';
export const REFERENCES_LIST_RESET_SELECTED =
  'REFERENCES--LIST--RESET_SELECTED';
export const REFERENCES_LIST_REQUEST_REFRESH =
  'REFERENCES--LIST--REQUEST_REFRESH';
export const REFERENCE_LOADING = 'REFERENCES--REFERENCE-LOADING';
export const REFERENCE_LOADED = 'REFERENCES--REFERENCE-LOADED';
export const REFERENCE_CACHED = 'REFERENCES--REFERENCE-CACHED';
export const REFERENCE_ERROR = 'REFERENCES--REFERENCE-ERROR';
export const REFERENCE_DELETED = 'REFERENCES--REFERENCE-DELETED';

export function setPerPage(perPage: number = 15) {
  return async (dispatch: Dispatch, getState: GetState) => {
    const {
      references: {
        referencesList: { per_page: oldPerPage }
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
        references: {
          referencesList: { refreshAll, refreshPages }
        }
      } = getState();
      if (refreshAll) dispatch(referencesListResetAll());
      if (refreshPages.length > 0)
        dispatch(referencesListResetSelected(refreshPages));
      dispatch(referencesListLoading());
      const {
        references: {
          referencesList: {
            pages,
            state: { per_page, sorting }
          }
        }
      } = getState();
      if (!has(pages, page)) {
        const references = await Api.References.fetch(per_page, sorting, page);
        dispatch(referencesListLoaded(references));
      } else {
        dispatch(referencesListCached(page));
      }
    } catch (e) {
      dispatch(referencesListError());
      dispatch(
        pushNotificationSimple(`An error occurred: ${e.message}!`, 'error')
      );
    }
  };
}

export function refreshPage(page: ?number = null) {
  return (dispatch: Dispatch) => {
    dispatch(referencesListRequestRefresh(!page ? null : [page]));
    if (page) dispatch(requestPage(page));
  };
}

export function requestReference(id: number, force: boolean = false) {
  return async (dispatch: Dispatch, getState: GetState) => {
    try {
      dispatch(referenceLoading());
      const {
        references: {
          references: { items }
        }
      } = getState();
      if (!has(items, id) || force) {
        const reference = await Api.References.fetchById(id);
        dispatch(referenceLoaded(reference));
      } else {
        dispatch(referenceCached());
      }
    } catch (e) {
      dispatch(referenceError());
      dispatch(
        pushNotificationSimple(`An error occurred: ${e.message}!`, 'error')
      );
    }
  };
}

export function deleteReference(id: number, page: ?number = null) {
  return async (dispatch: Dispatch) => {
    try {
      await Api.References.delete(id);
      dispatch(pushNotificationSimple('Reference deleted.'));
      dispatch(referenceDeleted(id));
      if (page) dispatch(refreshPage(page));
    } catch (e) {
      dispatch(
        pushNotificationSimple(`An error occurred: ${e.message}!`, 'error')
      );
    }
  };
}

function internalSetPerPage(perPage: number): Action {
  return {
    type: REFERENCES_LIST_SET_PER_PAGE,
    payload: {
      per_page: perPage
    }
  };
}

function internalSetSorting(payload: SortingSpec): Action {
  return {
    type: REFERENCES_LIST_SET_SORTING,
    payload
  };
}

export function referencesListLoading(): Action {
  return {
    type: REFERENCES_LIST_LOADING,
    payload: {}
  };
}

export function referencesListLoaded(payload: ReferencesCollection): Action {
  return {
    type: REFERENCES_LIST_LOADED,
    payload
  };
}

export function referencesListCached(page: number): Action {
  return {
    type: REFERENCES_LIST_CACHED,
    payload: {
      page
    }
  };
}

export function referencesListError(): Action {
  return {
    type: REFERENCES_LIST_ERROR,
    payload: {}
  };
}

export function referencesListResetAll(): Action {
  return {
    type: REFERENCES_LIST_RESET_ALL,
    payload: {}
  };
}

export function referencesListResetSelected(pages: number[]): Action {
  return {
    type: REFERENCES_LIST_RESET_SELECTED,
    payload: { pages }
  };
}

export function referencesListRequestRefresh(
  pages: ?(number[]) = null
): Action {
  return {
    type: REFERENCES_LIST_REQUEST_REFRESH,
    payload: {
      all: pages !== null,
      pages
    }
  };
}

export function referenceLoading(): Action {
  return {
    type: REFERENCE_LOADING,
    payload: {}
  };
}

export function referenceDeleted(payload: number): Action {
  return {
    type: REFERENCE_DELETED,
    payload
  };
}

export function referenceLoaded(payload: Reference): Action {
  return {
    type: REFERENCE_LOADED,
    payload
  };
}

export function referenceCached(): Action {
  return {
    type: REFERENCE_CACHED,
    payload: {}
  };
}

export function referenceError(): Action {
  return {
    type: REFERENCE_ERROR,
    payload: {}
  };
}
