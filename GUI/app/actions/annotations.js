/* eslint-disable camelcase */
// @flow
import { has } from 'lodash';
import type { Action, Dispatch, GetState } from '../reducers/types';
import type { Annotation, AnnotationsCollection } from '../types/annotations';
import * as Api from '../api';
import { pushNotificationSimple } from './notifications';
import type { SortingSpec } from '../types/common';

export const ANNOTATIONS_LIST_LOADING = 'ANNOTATIONS--LIST--LOADING';
export const ANNOTATIONS_LIST_LOADED = 'ANNOTATIONS--LIST--LOADED';
export const ANNOTATIONS_LIST_CACHED = 'ANNOTATIONS--LIST--CACHED';
export const ANNOTATIONS_LIST_ERROR = 'ANNOTATIONS--LIST--ERROR';
export const ANNOTATIONS_LIST_SET_PER_PAGE = 'ANNOTATIONS--LIST--SET_PER_PAGE';
export const ANNOTATIONS_LIST_SET_SORTING = 'ANNOTATIONS--LIST--SET_SORTING';
export const ANNOTATIONS_LIST_RESET_ALL = 'ANNOTATIONS--LIST--RESET_ALL';
export const ANNOTATIONS_LIST_RESET_SELECTED =
  'ANNOTATIONS--LIST--RESET_SELECTED';
export const ANNOTATIONS_LIST_REQUEST_REFRESH =
  'ANNOTATIONS--LIST--REQUEST_REFRESH';
export const ANNOTATION_LOADING = 'ANNOTATIONS--ANNOTATION-LOADING';
export const ANNOTATION_LOADED = 'ANNOTATIONS--ANNOTATION-LOADED';
export const ANNOTATION_CACHED = 'ANNOTATIONS--ANNOTATION-CACHED';
export const ANNOTATION_ERROR = 'ANNOTATIONS--ANNOTATION-ERROR';
export const ANNOTATION_DELETED = 'ANNOTATIONS--ANNOTATION-DELETED';

export function setPerPage(perPage: number = 15) {
  return async (dispatch: Dispatch, getState: GetState) => {
    const {
      annotations: {
        annotationsList: { per_page: oldPerPage }
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
        annotations: {
          annotationsList: { refreshAll, refreshPages }
        }
      } = getState();
      if (refreshAll) dispatch(annotationsListResetAll());
      if (refreshPages.length > 0)
        dispatch(annotationsListResetSelected(refreshPages));
      dispatch(annotationsListLoading());
      const {
        annotations: {
          annotationsList: {
            pages,
            state: { per_page, sorting }
          }
        }
      } = getState();
      if (!has(pages, page)) {
        const annotations = await Api.Annotations.fetch(
          per_page,
          sorting,
          page
        );
        dispatch(annotationsListLoaded(annotations));
      } else {
        dispatch(annotationsListCached(page));
      }
    } catch (e) {
      dispatch(annotationsListError());
      dispatch(
        pushNotificationSimple(`An error occurred: ${e.message}!`, 'error')
      );
    }
  };
}

export function refreshPage(page: ?number = null) {
  return (dispatch: Dispatch) => {
    dispatch(annotationsListRequestRefresh(!page ? null : [page]));
    if (page) dispatch(requestPage(page));
  };
}

export function requestAnnotation(id: number, force: boolean = false) {
  return async (dispatch: Dispatch, getState: GetState) => {
    try {
      dispatch(annotationLoading());
      const {
        annotations: {
          annotations: { items }
        }
      } = getState();
      if (!has(items, id) || force) {
        const annotation = await Api.Annotations.fetchById(id);
        dispatch(annotationLoaded(annotation));
      } else {
        dispatch(annotationCached());
      }
    } catch (e) {
      dispatch(annotationError());
      dispatch(
        pushNotificationSimple(`An error occurred: ${e.message}!`, 'error')
      );
    }
  };
}

export function deleteAnnotation(id: number, page: ?number = null) {
  return async (dispatch: Dispatch) => {
    try {
      await Api.Annotations.delete(id);
      dispatch(pushNotificationSimple('Annotation deleted.'));
      dispatch(annotationDeleted(id));
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
    type: ANNOTATIONS_LIST_SET_PER_PAGE,
    payload: {
      per_page: perPage
    }
  };
}

function internalSetSorting(payload: SortingSpec): Action {
  return {
    type: ANNOTATIONS_LIST_SET_SORTING,
    payload
  };
}

export function annotationsListLoading(): Action {
  return {
    type: ANNOTATIONS_LIST_LOADING,
    payload: {}
  };
}

export function annotationsListLoaded(payload: AnnotationsCollection): Action {
  return {
    type: ANNOTATIONS_LIST_LOADED,
    payload
  };
}

export function annotationsListCached(page: number): Action {
  return {
    type: ANNOTATIONS_LIST_CACHED,
    payload: {
      page
    }
  };
}

export function annotationsListError(): Action {
  return {
    type: ANNOTATIONS_LIST_ERROR,
    payload: {}
  };
}

export function annotationsListResetAll(): Action {
  return {
    type: ANNOTATIONS_LIST_RESET_ALL,
    payload: {}
  };
}

export function annotationsListResetSelected(pages: number[]): Action {
  return {
    type: ANNOTATIONS_LIST_RESET_SELECTED,
    payload: { pages }
  };
}

export function annotationsListRequestRefresh(
  pages: ?(number[]) = null
): Action {
  return {
    type: ANNOTATIONS_LIST_REQUEST_REFRESH,
    payload: {
      all: pages !== null,
      pages
    }
  };
}

export function annotationLoading(): Action {
  return {
    type: ANNOTATION_LOADING,
    payload: {}
  };
}

export function annotationDeleted(payload: number): Action {
  return {
    type: ANNOTATION_DELETED,
    payload
  };
}

export function annotationLoaded(payload: Annotation): Action {
  return {
    type: ANNOTATION_LOADED,
    payload
  };
}

export function annotationCached(): Action {
  return {
    type: ANNOTATION_CACHED,
    payload: {}
  };
}

export function annotationError(): Action {
  return {
    type: ANNOTATION_ERROR,
    payload: {}
  };
}
