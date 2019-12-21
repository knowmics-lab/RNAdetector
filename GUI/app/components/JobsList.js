// @flow
import React, { Component } from 'react';
import has from 'lodash/has';
import { Typography, CircularProgress } from '@material-ui/core';
import Table from './UI/Table';
import type { JobsCollectionItem, JobsListType } from '../types/jobs';
import Snackbar from './UI/Snackbar';
import type { StatePaginationType } from '../types/common';

type Props = {
  requestPage: number => void,
  setPerPage: number => void,
  jobsListResetLoading: () => void,
  refreshAll: boolean,
  refreshPages: number[],
  state: StatePaginationType,
  pages: { +[number]: JobsCollectionItem[] }
} & JobsListType;

export default class JobsList extends Component<Props> {
  props: Props;

  handleClose = () => {
    const { jobsListResetLoading } = this.props;
    jobsListResetLoading();
  };

  render() {
    const {
      refreshAll,
      refreshPages,
      state,
      pages,
      requestPage,
      setPerPage
    } = this.props;
    const cp = state.current_page;
    const needsLoading =
      !state.isError &&
      !state.isFetching &&
      (!cp || refreshAll || refreshPages.includes(cp) || !has(pages, cp));
    if (needsLoading) {
      requestPage(1);
    }
    return (
      <div>
        <Typography paragraph>Hello World</Typography>
        {needsLoading ? (
          <div className="text-center">
            <CircularProgress />
          </div>
        ) : (
          !state.isError && (
            <Table
              columns={[
                {
                  id: 'type',
                  label: 'Type',
                  minWidth: 100
                },
                {
                  id: 'status',
                  label: 'Status',
                  minWidth: 100
                },
                {
                  id: 'created_at',
                  label: 'Created at',
                  minWidth: 150
                }
              ]}
              currentPage={state.current_page || 1}
              data={pages[state.current_page] || null}
              keyField="id"
              loading={state.isFetching}
              rowsPerPage={state.per_page}
              totalRows={state.total || 1}
              onPageChange={page => requestPage(page)}
              onRowsPerPageChange={perPage => setPerPage(perPage)}
            />
          )
        )}
        <Snackbar
          message={`An error occurred: ${state.errorMessage}!`}
          isOpen={state.isError}
          setClosed={this.handleClose}
          variant="error"
        />
      </div>
    );
  }
}
