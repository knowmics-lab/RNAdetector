/* eslint-disable no-case-declarations */
// @flow
import * as AnnotationsActions from '../actions/annotations';
import { Action } from './types';
import * as Api from '../api';
import type {
  AnnotationsListType,
  AnnotationsStateType,
  AnnotationsCollection,
  LoadedAnnotations
} from '../types/annotations';
import type { SortingSpec } from '../types/common';

const initState = (
  perPage: number = 15,
  sorting: SortingSpec = { created_at: 'desc' }
): AnnotationsStateType => ({
  annotationsList: {
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
  annotations: {
    fetching: false,
    items: {}
  }
});

function resetList(
  state: AnnotationsListType,
  pages: number[]
): AnnotationsListType {
  return {
    ...state,
    refreshPages: [],
    pages: Api.Utils.filterByKey(
      state.pages,
      key => !pages.includes(parseInt(key, 10))
    )
  };
}

function changeListState(state, newState) {
  return {
    ...state,
    state: {
      ...state.state,
      ...newState
    }
  };
}

function addLoadedPayload(
  state: AnnotationsListType,
  payload: AnnotationsCollection
): AnnotationsListType {
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

function handleList(
  state: AnnotationsListType,
  action: Action
): AnnotationsListType {
  switch (action.type) {
    case AnnotationsActions.ANNOTATIONS_LIST_RESET_ALL:
      const newState = initState(state.state.per_page, state.state.sorting);
      return newState.annotationsList;
    case AnnotationsActions.ANNOTATIONS_LIST_RESET_SELECTED:
      return resetList(state, action.payload.pages);
    case AnnotationsActions.ANNOTATIONS_LIST_ERROR:
      return changeListState(state, {
        current_page: state.state.current_page || 0,
        last_page: state.state.last_page || 0,
        total: state.state.total || 0,
        fetching: false
      });
    case AnnotationsActions.ANNOTATIONS_LIST_SET_PER_PAGE:
      return {
        ...changeListState(state, {
          per_page: action.payload.per_page
        }),
        refreshAll: true,
        refreshPages: []
      };
    case AnnotationsActions.ANNOTATIONS_LIST_SET_SORTING:
      return {
        ...changeListState(state, { sorting: action.payload }),
        refreshAll: true,
        refreshPages: []
      };
    case AnnotationsActions.ANNOTATIONS_LIST_LOADING:
      return changeListState(state, {
        fetching: true
      });
    case AnnotationsActions.ANNOTATIONS_LIST_REQUEST_REFRESH:
      return {
        ...state,
        refreshAll: state.refreshAll || action.payload.all,
        refreshPages: [...state.refreshPages, ...(action.payload.pages || [])]
      };
    case AnnotationsActions.ANNOTATIONS_LIST_LOADED:
      return addLoadedPayload(state, action.payload);
    case AnnotationsActions.ANNOTATIONS_LIST_CACHED:
      return changeListState(state, {
        current_page: action.payload.page,
        fetching: false
      });
    default:
      return state;
  }
}

function handleFetching(state, fetching) {
  return {
    ...state,
    fetching
  };
}

function handleOne(
  state: LoadedAnnotations,
  action: Action
): LoadedAnnotations {
  switch (action.type) {
    case AnnotationsActions.ANNOTATION_ERROR:
    case AnnotationsActions.ANNOTATION_LOADING:
      return handleFetching(state, true);
    case AnnotationsActions.ANNOTATION_LOADED:
      return {
        fetching: false,
        items: {
          ...state.items,
          [action.payload.id]: action.payload
        }
      };
    case AnnotationsActions.ANNOTATION_CACHED:
      return handleFetching(state, false);
    case AnnotationsActions.ANNOTATION_DELETED:
      return {
        fetching: state.fetching,
        items: Api.Utils.filterByKey(state.items, k => k !== action.payload)
      };
    default:
      return state;
  }
}

export default function annotations(
  state: AnnotationsStateType,
  action: Action
): AnnotationsStateType {
  const oldState = state || initState();
  const newState = {
    annotationsList: handleList(oldState.annotationsList, action),
    annotations: handleOne(oldState.annotations, action)
  };
  if (
    newState.annotationsList === oldState.annotationsList &&
    newState.annotations === oldState.annotations
  ) {
    return oldState;
  }
  return newState;
}
