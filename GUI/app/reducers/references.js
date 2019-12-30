/* eslint-disable no-case-declarations */
// @flow
import * as ReferencesActions from '../actions/references';
import { Action } from './types';
import * as Api from '../api';
import type {
  ReferencesListType,
  ReferencesStateType,
  ReferencesCollection,
  LoadedReferences
} from '../types/references';
import type { SortingSpec } from '../types/common';

const initState = (
  perPage: number = 15,
  sorting: SortingSpec = { created_at: 'desc' }
): ReferencesStateType => ({
  referencesList: {
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
  references: {
    fetching: false,
    items: {}
  }
});

function resetList(
  state: ReferencesListType,
  pages: number[]
): ReferencesListType {
  return {
    ...state,
    refreshPages: [],
    pages: Api.Utils.filterByKey(state.pages, key => !pages.includes(key))
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
  state: ReferencesListType,
  payload: ReferencesCollection
): ReferencesListType {
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
  state: ReferencesListType,
  action: Action
): ReferencesListType {
  switch (action.type) {
    case ReferencesActions.REFERENCES_LIST_RESET_ALL:
      const newState = initState(state.state.per_page, state.state.sorting);
      return newState.referencesList;
    case ReferencesActions.REFERENCES_LIST_RESET_SELECTED:
      return resetList(state, action.payload.pages);
    case ReferencesActions.REFERENCES_LIST_ERROR:
      return changeListState(state, {
        current_page: state.state.current_page || 0,
        last_page: state.state.last_page || 0,
        total: state.state.total || 0,
        fetching: false
      });
    case ReferencesActions.REFERENCES_LIST_SET_PER_PAGE:
      return {
        ...changeListState(state, {
          per_page: action.payload.per_page
        }),
        refreshAll: true,
        refreshPages: []
      };
    case ReferencesActions.REFERENCES_LIST_SET_SORTING:
      return {
        ...changeListState(state, { sorting: action.payload }),
        refreshAll: true,
        refreshPages: []
      };
    case ReferencesActions.REFERENCES_LIST_LOADING:
      return changeListState(state, {
        fetching: true
      });
    case ReferencesActions.REFERENCES_LIST_REQUEST_REFRESH:
      return {
        ...state,
        refreshAll: state.refreshAll || action.payload.all,
        refreshPages: [...state.refreshPages, ...(action.payload.pages || [])]
      };
    case ReferencesActions.REFERENCES_LIST_LOADED:
      return addLoadedPayload(state, action.payload);
    case ReferencesActions.REFERENCES_LIST_CACHED:
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

function handleOne(state: LoadedReferences, action: Action): LoadedReferences {
  switch (action.type) {
    case ReferencesActions.REFERENCE_ERROR:
    case ReferencesActions.REFERENCE_LOADING:
      return handleFetching(state, true);
    case ReferencesActions.REFERENCE_LOADED:
      return {
        fetching: false,
        items: {
          ...state.items,
          [action.payload.id]: action.payload
        }
      };
    case ReferencesActions.REFERENCE_CACHED:
      return handleFetching(state, false);
    case ReferencesActions.REFERENCE_DELETED:
      return {
        fetching: state.fetching,
        items: Api.Utils.filterByKey(state.items, k => k !== action.payload)
      };
    default:
      return state;
  }
}

export default function references(
  state: ReferencesStateType,
  action: Action
): ReferencesStateType {
  const oldState = state || initState();
  const newState = {
    referencesList: handleList(oldState.referencesList, action),
    references: handleOne(oldState.references, action)
  };
  if (
    newState.referencesList === oldState.referencesList &&
    newState.references === oldState.references
  ) {
    return oldState;
  }
  return newState;
}
