/* eslint-disable react/no-unused-prop-types,no-unused-vars,camelcase */
// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import Paper from '@material-ui/core/Paper';
import Box from '@material-ui/core/Box';
import FormGroup from '@material-ui/core/FormGroup';
import Grid from '@material-ui/core/Grid';
import CircularProgress from '@material-ui/core/CircularProgress';
import { Formik, Form, FieldArray } from 'formik';
import * as Yup from 'yup';
import Backdrop from '@material-ui/core/Backdrop';
import { InputLabel } from '@material-ui/core';
import Button from '@material-ui/core/Button';
import Icon from '@material-ui/core/Icon';
import IconButton from '@material-ui/core/IconButton';
import * as Api from '../../api';
import * as DiffExpConsts from '../../constants/diff-exp-consts';
import * as DiffExpDefaults from '../../constants/diff-exp-defaults';
import { JOBS } from '../../constants/routes';
import SelectField from '../Form/SelectField';
import TextField from '../Form/TextField';
import Wizard from '../UI/Wizard';
import SwitchField from '../Form/SwitchField';
import type { Job } from '../../types/jobs';
import TableField from '../Form/TableField';
import type { PushNotificationFunction } from '../../types/notifications';
import { SubmitButton } from '../UI/Button';
import FileSelector from '../UI/FileSelector';

type Props = {
  refreshJobs: () => void,
  redirect: mixed => void,
  pushNotification: PushNotificationFunction,
  classes: {
    root: *,
    formControl: *,
    buttonWrapper: *,
    buttonProgress: *,
    backButton: *,
    instructions: *,
    backdrop: *
  }
};

const style = theme => ({
  root: {
    padding: theme.spacing(3, 2)
  },
  formControl: {
    margin: theme.spacing(1),
    minWidth: 120
  },
  backButton: {
    marginRight: theme.spacing(1)
  },
  instructions: {
    marginTop: theme.spacing(1),
    marginBottom: theme.spacing(1)
  },
  backdrop: {
    zIndex: theme.zIndex.drawer + 1,
    color: '#fff'
  }
});

const SUPPORTED_ANALYSIS = {
  long_rna_job_type: 'Long RNA',
  small_rna_job_type: 'Small RNA',
  circ_rna_job_type: 'Circ RNA'
};

type State = {
  isLoading: boolean,
  isSaving: boolean,
  validationErrors: *,
  selectedJob: ?Job,
  variables: { [string]: string },
  variablesMeta: { [string]: string[] }
};

class DiffExpr extends React.Component<Props, State> {
  props: Props;

  cachedContrasts: ?{
    key: ?string,
    contrasts: ?{ [string]: string }
  };

  constructor(props) {
    super(props);
    this.cachedContrasts = null;
    this.state = {
      isLoading: false,
      isSaving: false,
      validationErrors: {},
      selectedJob: null,
      variables: {},
      variablesMeta: {}
    };
  }

  loadJob = id => {
    if (!id) {
      this.cachedContrasts = null;
      this.setState({
        selectedJob: null,
        variables: {},
        variablesMeta: {}
      });
    } else {
      this.setState({
        isLoading: true
      });
      const { pushNotification } = this.props;
      Api.Jobs.fetchJobById(id)
        // eslint-disable-next-line promise/always-return
        .then(selectedJob => {
          this.cachedContrasts = null;
          this.setState({
            isLoading: false,
            selectedJob,
            variables: this.getVariables(selectedJob),
            variablesMeta: this.getMeta(selectedJob)
          });
        })
        .catch(e => {
          pushNotification(`An error occurred: ${e.message}!`, 'error');
          this.setState({
            isLoading: false
          });
        });
    }
  };

  // eslint-disable-next-line class-methods-use-this
  getMeta(selectedJob: Job): { [string]: string[] } {
    if (
      !selectedJob ||
      !selectedJob.output ||
      !Array.isArray(selectedJob.output.metadata)
    ) {
      return {};
    }
    return Object.fromEntries(
      selectedJob.output.metadata
        .filter(f => f.type === 'string')
        // $FlowFixMe
        .map(f => [f.field, f.content])
    );
  }

  // eslint-disable-next-line class-methods-use-this
  getVariables(selectedJob: Job): { [string]: string } {
    if (
      !selectedJob ||
      !selectedJob.output ||
      !Array.isArray(selectedJob.output.metadata)
    ) {
      return {};
    }
    return Object.fromEntries(
      selectedJob.output.metadata
        // $FlowFixMe
        .map(f => (f.type === 'string' ? f.field : null))
        .filter(v => v !== null)
        // $FlowFixMe
        .map(v => [v, v])
    );
  }

  recursiveConditionBuilder(variables: string[], prefix: string) {
    const { variablesMeta } = this.state;
    const [first, ...rest] = variables;
    let results: [string, string][] = [];
    variablesMeta[first].forEach(v => {
      const elem = `${prefix}${v}`;
      if (rest.length === 0) {
        results.push([elem, elem]);
      } else {
        results = [
          ...results,
          ...this.recursiveConditionBuilder(rest, `${elem}_`)
        ];
      }
    });
    return results;
  }

  getContrasts(conditionVariables: ?(string[])): { [string]: string } {
    if (!conditionVariables) return {};
    if (conditionVariables.length <= 0) return {};
    const { selectedJob } = this.state;
    if (!selectedJob) return {};
    const key = JSON.stringify(conditionVariables);
    if (
      this.cachedContrasts &&
      this.cachedContrasts.key === key &&
      this.cachedContrasts.contrasts
    ) {
      return this.cachedContrasts.contrasts;
    }
    const contrasts = Object.fromEntries(
      this.recursiveConditionBuilder(conditionVariables, '')
    );
    this.cachedContrasts = {
      key,
      contrasts
    };
    return contrasts;
  }

  getValidationSchema = () =>
    Yup.object().shape({
      /* code: Yup.string()
        .required()
        .matches(/^[A-Za-z0-9\-_]+$/, {
          message: 'The field must contain only letters, numbers, and dashes.'
        }),
      name: Yup.string().required(),
      de_novo: Yup.boolean().notRequired(),
      analysis: Yup.string()
        .required()
        .oneOf(Object.keys(SUPPORTED_ANALYSIS)),
      jobs: Yup.array()
        .required()
        .min(1)
        .of(Yup.number()) */
    });

  getSteps = () => [
    'Choose name and type',
    'Select Variables and Contrasts',
    'Common Parameters',
    'Normalization Parameters',
    'Statistics Parameters',
    'Start Analysis'
  ];

  getStep0 = () => {
    const { classes, pushNotification } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can choose code and name for the analysis, and select a
          sample group that will be used for the analysis. Sample code is used
          for further analysis. It identifies the sample during the analysis
          therefore it should be a string without any spaces (only letters,
          numbers, and dashes).
        </Typography>
        <TextField label="Code" name="code" required />
        <TextField label="Analysis Name" name="name" required />
        <TableField
          name="source_sample_group"
          required
          single
          getData={() => Api.Jobs.fetchAllByType('samples_group_job_type')}
          onError={e =>
            pushNotification(`An error occurred: ${e.message}!`, 'error')
          }
          onChange={id => this.loadJob(id)}
          label="Select one sample group"
          columns={[
            {
              dataField: 'sample_code',
              label: 'Code'
            },
            {
              dataField: 'name',
              label: 'Name'
            },
            {
              dataField: 'created_at_diff',
              sortingField: 'created_at',
              label: 'Created at'
            }
          ]}
        />
      </>
    );
  };

  addHandle = helpers => e => {
    e.preventDefault();
    helpers.push({ case: '', control: '' });
  };

  removeHandle = (helpers, i) => e => {
    e.preventDefault();
    helpers.remove(i);
  };

  makeSampleForm = (i, single, values, helpers) => {
    const { classes } = this.props;
    return (
      <FormGroup row className={classes.formControl} key={`sample-${i}`}>
        <Grid container justify="space-around" alignItems="center" spacing={3}>
          <Grid item xs>
            <SelectField
              label="Case"
              name={`contrasts.${i}.case`}
              options={values}
              addEmpty
            />
          </Grid>
          <Grid item xs>
            <SelectField
              label="Case"
              name={`contrasts.${i}.control`}
              options={values}
              addEmpty
            />
          </Grid>
          {!single && (
            <Grid item xs={1}>
              <IconButton onClick={this.removeHandle(helpers, i)}>
                <Icon className="fas fa-trash" />
              </IconButton>
            </Grid>
          )}
        </Grid>
      </FormGroup>
    );
  };

  getStep1 = values => {
    const { classes } = this.props;
    const { variables } = this.state;
    const { condition_variables, contrasts } = values;
    const availableValues = this.getContrasts(condition_variables);
    const emptyValues =
      availableValues && Object.keys(availableValues).length === 0;
    const single = !contrasts || contrasts.length <= 1;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can choose which variables will be used to define the
          contrasts. After selecting the variables, you will be able to add one
          or more contrasts each with a case and a control.
        </Typography>
        <SelectField
          label="Select condition variables:"
          name="condition_variables"
          options={variables}
          required
          multiple
        />
        {!emptyValues && (
          <FieldArray
            name="contrasts"
            render={helpers => (
              <>
                {contrasts &&
                  contrasts.map((s, i) =>
                    this.makeSampleForm(i, single, availableValues, helpers)
                  )}
                <FormGroup row className={classes.formControl}>
                  <Grid
                    container
                    direction="row-reverse"
                    alignItems="center"
                    spacing={3}
                  >
                    <Grid item xs="auto">
                      <Button
                        variant="outlined"
                        onClick={this.addHandle(helpers)}
                      >
                        <Icon className="fas fa-plus" /> Add a contrast
                      </Button>
                    </Grid>
                  </Grid>
                </FormGroup>
              </>
            )}
          />
        )}
      </>
    );
  };

  getStep2 = () => {
    const { classes } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can upload an optional sample description table (in TSV
          format) that can be used for samples annotation and differential
          expression analysis. Finally, to proceed with the creation of a new
          sample group, click on the &quot;Save&quot; button.
        </Typography>
        <TextField
          label="p-value Threshold"
          name="parameters.pcut"
          required
          type="number"
          helperText="A p-value cutoff for exporting differentially genes."
        />
        <TextField
          label="Log offset"
          name="parameters.log_offset"
          required
          type="number"
          helperText="An offset to be added to values during logarithmic transformations in order to avoid Infinity."
        />
        <SelectField
          label="When filtering should be applied?"
          name="parameters.when_apply_filter"
          options={DiffExpConsts.when_apply_filter}
        />
        <SelectField
          label="Which method should be used for p-value adjustment?"
          name="parameters.adjust_method"
          options={DiffExpConsts.adjust_method}
          required
        />
        <SelectField
          label="Which method should be used for p-value meta-analysis?"
          name="parameters.meta_p_method"
          options={DiffExpConsts.meta_p_method}
          required
        />
        <SelectField
          label="Which formats should be used to export figures?"
          name="parameters.fig_formats"
          options={DiffExpConsts.fig_formats}
          multiple
          required
        />
        <TextField
          label="Number of threads"
          name="parameters.num_cores"
          type="number"
          helperText={`Do not select more than ${Math.floor(
            Api.Utils.cpuCount() / 3
          )} cores if you wish to submit multiple analysis.`}
          required
        />
      </>
    );
  };

  getStep3 = () => {};

  getStep4 = () => {};

  getStep5 = () => {};

  getSubmitButton = () => {
    const { isSaving } = this.state;
    return <SubmitButton isSaving={isSaving} />;
  };

  setSaving = (isSaving, validationErrors = {}) => {
    this.setState({
      isSaving,
      validationErrors
    });
  };

  formSubmit = async values => {
    console.log(values);
    /* const { pushNotification, redirect, refreshJobs } = this.props;
    this.setSaving(true);
    try {
      const groupJob = await this.createGroup(values);
      if (groupJob) {
        await Api.Jobs.submitJob(groupJob.id);
        refreshJobs();
        redirect(JOBS);
      }
    } catch (e) {
      pushNotification(`An error occurred: ${e.message}`, 'error');
      this.setSaving(false);
    } */
  };

  render() {
    const { classes } = this.props;
    const { validationErrors, isLoading } = this.state;
    const steps = this.getSteps();
    return (
      <>
        <Box>
          <Paper className={classes.root}>
            <Typography variant="h5" component="h3">
              Differential Expression Analysis
            </Typography>
            <Formik
              initialValues={DiffExpDefaults.defaults}
              initialErrors={validationErrors}
              validationSchema={this.getValidationSchema()}
              onSubmit={v => {
                this.formSubmit(v).catch(() => false);
              }}
            >
              {({ values }) => (
                <Form>
                  <Wizard steps={steps} submitButton={this.getSubmitButton}>
                    <div>{this.getStep0()}</div>
                    <div>{this.getStep1(values)}</div>
                    <div>{this.getStep2()}</div>
                    <div>{this.getStep3()}</div>
                    <div>{this.getStep4()}</div>
                    <div>{this.getStep5()}</div>
                  </Wizard>
                </Form>
              )}
            </Formik>
          </Paper>
        </Box>
        <Backdrop className={classes.backdrop} open={isLoading}>
          <CircularProgress color="inherit" />
        </Backdrop>
      </>
    );
  }
}

export default withStyles(style)(DiffExpr);
