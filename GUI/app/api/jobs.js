// @flow

export type JobType = {
  id: number,
  type: string,
  status: 'ready' | 'queued' | 'processing' | 'completed' | 'failed',
  parameters: { [string]: string | number },
  output: { [string]: string | number },
  log: string,
  created_at: string,
  updated_at: string,
  owner: *,
  links: {
    self: string,
    owner: string,
    upload: string,
    submit: string
  }
};

export type JobsCollection = {};

export default {

};
