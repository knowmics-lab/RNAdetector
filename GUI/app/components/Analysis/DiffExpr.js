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
import { Formik, Form, FieldArray, Field } from 'formik';
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
import type {
  ContrastType,
  DiffExpParameters,
  SampleTypes
} from '../../types/analysis';

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
    instructionsSmall: *,
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
  instructionsSmall: {
    margin: theme.spacing(1),
    fontSize: '0.8rem'
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

  hasTranscripts(): boolean {
    const { selectedJob } = this.state;
    if (!selectedJob) return false;
    if (!selectedJob.output) return false;
    return !!selectedJob.output.harmonizedTranscriptsFile;
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

  sortConditionVariables(selectedVars: string[]): string[] {
    const { variables } = this.state;
    return Object.keys(variables).filter(v => selectedVars.includes(v));
  }

  getContrasts(conditionVariables: ?(string[])): { [string]: string } {
    if (!conditionVariables) return {};
    if (conditionVariables.length <= 0) return {};
    const { selectedJob } = this.state;
    if (!selectedJob) return {};
    const sortedVariables = this.sortConditionVariables(conditionVariables);
    const key = JSON.stringify(sortedVariables);
    if (
      this.cachedContrasts &&
      this.cachedContrasts.key === key &&
      this.cachedContrasts.contrasts
    ) {
      return this.cachedContrasts.contrasts;
    }
    const contrasts = Object.fromEntries(
      this.recursiveConditionBuilder(sortedVariables, '')
    );
    this.cachedContrasts = {
      key,
      contrasts
    };
    return contrasts;
  }

  getValidationSchema = () =>
    Yup.object().shape({
      code: Yup.string()
        .required()
        .matches(/^[A-Za-z0-9\-_]+$/, {
          message: 'The field must contain only letters, numbers, and dashes.'
        }),
      name: Yup.string().required(),
      source_sample_group: Yup.number().required(),
      sample_type: Yup.string().oneOf(Object.keys(DiffExpConsts.sample_type)),
      condition_variables: Yup.array()
        .of(Yup.string())
        .min(1),
      contrasts: Yup.array()
        .of(
          Yup.object().shape({
            case: Yup.string().required('You must select one case'),
            control: Yup.string().required('You must select one control')
          })
        )
        .min(1),
      parameters: Yup.object().shape({
        pcut: Yup.number()
          .positive()
          .min(0)
          .max(1),
        log_offset: Yup.number().positive(),
        when_apply_filter: Yup.string().oneOf(
          Object.keys(DiffExpConsts.when_apply_filter)
        ),
        norm: Yup.string().oneOf(Object.keys(DiffExpConsts.norm)),
        norm_args: Yup.object().shape({
          method: Yup.string().oneOf(
            Object.keys(DiffExpConsts.norm_args_method)
          ),
          locfunc: Yup.string().oneOf(
            Object.keys(DiffExpConsts.norm_args_locfunc)
          )
        }),
        stats: Yup.array().of(
          Yup.string().oneOf(Object.keys(DiffExpConsts.stats))
        ),
        stats_args: Yup.object().shape({
          deseq: Yup.object().shape({
            fitType: Yup.string().oneOf(
              Object.keys(DiffExpConsts.stats_args.deseq.fitType)
            )
          }),
          edger: Yup.object().shape({
            main_method: Yup.string().oneOf(
              Object.keys(DiffExpConsts.stats_args.edger.main_method)
            ),
            rowsum_filter: Yup.number().integer(),
            trend: Yup.string().oneOf(
              Object.keys(DiffExpConsts.stats_args.edger.trend)
            ),
            tag_method: Yup.string().oneOf(
              Object.keys(DiffExpConsts.stats_args.edger.tag_method)
            ),
            glm_method: Yup.string().oneOf(
              Object.keys(DiffExpConsts.stats_args.edger.glm_method)
            ),
            trend_method: Yup.string().oneOf(
              Object.keys(DiffExpConsts.stats_args.edger.trend_method)
            )
          }),
          limma: Yup.object().shape({
            normalize_method: Yup.string().oneOf(
              Object.keys(DiffExpConsts.stats_args.limma.normalize_method)
            )
          })
        }),
        filters: Yup.object().shape({
          enabled: Yup.array().of(
            Yup.string().oneOf(Object.keys(DiffExpConsts.filters))
          ),
          length: Yup.object().shape({
            length: Yup.number().integer()
          }),
          avg_reads: Yup.object().shape({
            average_per_bp: Yup.number().integer(),
            quantile: Yup.number()
              .min(0)
              .max(1)
          }),
          expression: Yup.object().shape({
            median: Yup.boolean(),
            mean: Yup.boolean(),
            quantile: Yup.number()
              .min(0)
              .max(1),
            known: Yup.array().of(Yup.string())
          }),
          presence: Yup.object().shape({
            frac: Yup.number()
              .min(0)
              .max(1),
            min_count: Yup.number().integer(),
            per_condition: Yup.boolean()
          })
        }),
        adjust_method: Yup.string().oneOf(
          Object.keys(DiffExpConsts.adjust_method)
        ),
        meta_p_method: Yup.string().oneOf(
          Object.keys(DiffExpConsts.meta_p_method)
        ),
        fig_formats: Yup.array()
          .of(Yup.string().oneOf(Object.keys(DiffExpConsts.fig_formats)))
          .min(1),
        num_cores: Yup.number().integer()
      })
    });

  getSteps = () => [
    'Choose name and type',
    'Select Variables and Contrasts',
    'Common Parameters',
    'Normalization Parameters',
    'Filtering Parameters',
    'Statistics Parameters'
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
    const hasTranscripts = this.hasTranscripts();
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can select the common parameters for the analysis.
        </Typography>
        {!hasTranscripts && <Field type="hidden" name="sample_type" />}
        {hasTranscripts && (
          <SelectField
            label="Which type of analysis you wish to run?"
            name="sample_type"
            options={DiffExpConsts.sample_type}
          />
        )}
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

  getStep3 = values => {
    const { classes } = this.props;
    const {
      parameters: { norm }
    } = values;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can select the normalization method and its parameters.
        </Typography>
        <SelectField
          label="Which normalization algorithm should be used?"
          name="parameters.norm"
          options={DiffExpConsts.norm}
          required
        />
        {norm === 'edger' && (
          <SelectField
            label="Normalization method"
            name="parameters.norm_args.method"
            options={DiffExpConsts.norm_args_method}
            required
          />
        )}
        {norm === 'deseq' && (
          <SelectField
            label="Location function"
            name="parameters.norm_args.locfunc"
            options={DiffExpConsts.norm_args_locfunc}
            required
            helperText="A function to compute a location for a sample. For low counts, the shorth function may give better results."
          />
        )}
      </>
    );
  };

  getStep4 = values => {
    const { classes } = this.props;
    const {
      parameters: {
        filters: {
          enabled,
          avg_reads: { average_per_bp },
          presence: { frac, min_count, per_condition }
        }
      }
    } = values;
    const hasFilt = f => enabled.includes(f);
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can determine how reads will be filtered from the raw data.
        </Typography>
        <SelectField
          label="Which type of filters are enabled?"
          name="parameters.filters.enabled"
          options={DiffExpConsts.filters}
          multiple
          required
        />
        {hasFilt('length') && (
          <Box mx={2} mt={2}>
            <Typography variant="subtitle2">Length filter</Typography>
            <TextField
              label="Length"
              name="parameters.filters.length.length"
              required
              type="number"
              helperText="Genes/transcripts are accepted for further analysis if they are above the length in kb"
            />
          </Box>
        )}
        {hasFilt('reads') && (
          <Box mx={2} mt={2}>
            <Typography variant="subtitle2">Reads filter</Typography>
            <Typography className={classes.instructionsSmall}>
              A gene/transcript is accepted for further analysis if it has more
              average reads than the specified quantile of the average count
              distribution per {average_per_bp} base pairs
            </Typography>
            <TextField
              label="Number of base pairs for the average"
              name="parameters.filters.avg_reads.average_per_bp"
              required
              type="number"
            />
            <TextField
              label="Quantile"
              name="parameters.filters.avg_reads.quantile"
              required
              type="number"
            />
          </Box>
        )}
        {hasFilt('reads') && (
          <Box mx={2} mt={2}>
            <Typography variant="subtitle2">Expression filter</Typography>
            <Typography className={classes.instructionsSmall}>
              A filter based on the overall expression of a gene
            </Typography>
            <SwitchField
              label="Median (Genes below the median of the overall count distribution are not accepted)"
              name="parameters.filters.expression.median"
              required
            />
            <SwitchField
              label="Mean  (Genes below the mean of the overall count distribution are not accepted)"
              name="parameters.filters.expression.mean"
              required
            />
            <TextField
              label="Quantile"
              name="parameters.filters.expression.quantile"
              type="number"
              helperText="Genes below the specified quantile of the overall count distribution are not accepted (leave empty to disable)."
            />
          </Box>
        )}
        {hasFilt('presence') && (
          <Box mx={2} mt={2}>
            <Typography variant="subtitle2">Presence filter</Typography>
            <Typography className={classes.instructionsSmall}>
              A gene is considered for statistical testing if {frac * 100}% of
              samples {per_condition ? 'per condition' : ''} have more than{' '}
              {min_count} reads
            </Typography>
            <TextField
              label="Fraction of samples"
              name="parameters.filters.presence.frac"
              type="number"
              helperText="Percentage should be in the range [0,1]."
            />
            <TextField
              label="Minimum number of reads"
              name="parameters.filters.presence.min_count"
              type="number"
            />
            <SwitchField
              label="Check presence per condition?"
              name="parameters.filters.presence.per_condition"
              required
            />
          </Box>
        )}
      </>
    );
  };

  getStep5 = values => {
    const { classes } = this.props;
    const {
      parameters: {
        stats,
        stats_args: {
          edger: { main_method }
        }
      }
    } = values;
    const hasStat = s => stats.includes(s);
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can select one or more statistical methods, and their
          parameters, to use for analysis. After entering all the parameters,
          you can start the analysis by clicking on the &quot;Start
          Analysis&quot; button.
        </Typography>
        <SelectField
          label="Which statistical methods should be used?"
          name="parameters.stats"
          options={DiffExpConsts.stats}
          multiple
          required
        />
        {hasStat('limma') && (
          <Box mx={2} mt={2}>
            <Typography variant="subtitle2">Limma parameters</Typography>
            <SelectField
              label="Voom normalize method"
              name="parameters.stats_args.limma.normalize_method"
              options={DiffExpConsts.stats_args.limma.normalize_method}
              required
            />
          </Box>
        )}
        {hasStat('deseq') && (
          <Box mx={2} mt={2}>
            <Typography variant="subtitle2">DESeq2 parameters</Typography>
            <SelectField
              label="Dispersion fitting method"
              name="parameters.stats_args.deseq.fitType"
              options={DiffExpConsts.stats_args.deseq.fitType}
              required
              helperText="The type of fitting of dispersions to the mean intensity."
            />
          </Box>
        )}
        {hasStat('edger') && (
          <Box mx={2} mt={2}>
            <Typography variant="subtitle2">edgeR parameters</Typography>
            <SelectField
              label="Type of analysis"
              name="parameters.stats_args.edger.main_method"
              options={DiffExpConsts.stats_args.edger.main_method}
              required
            />
            {main_method === 'classic' && (
              <>
                <TextField
                  label="Row sum filter"
                  name="parameters.stats_args.edger.rowsum_filter"
                  required
                  type="number"
                  helperText="Genes with total count (across all samples) below this value will be filtered out before estimating the dispersion."
                />
                <SelectField
                  label="Dispersion trend estimation method"
                  name="parameters.stats_args.edger.trend"
                  options={DiffExpConsts.stats_args.edger.trend}
                  required
                />
                <SelectField
                  label="Posterior likelihood maximization method"
                  name="parameters.stats_args.edger.tag_method"
                  options={DiffExpConsts.stats_args.edger.tag_method}
                  required
                />
              </>
            )}
            {main_method === 'glm' && (
              <>
                <SelectField
                  label="Dispersion estimation method"
                  name="parameters.stats_args.edger.glm_method"
                  options={DiffExpConsts.stats_args.edger.glm_method}
                  required
                />
                <SelectField
                  label="Trended dispersion estimation method"
                  name="parameters.stats_args.edger.trend_method"
                  options={DiffExpConsts.stats_args.edger.trend_method}
                  required
                />
              </>
            )}
          </Box>
        )}
      </>
    );
  };

  getSubmitButton = () => {
    const { isSaving } = this.state;
    return <SubmitButton text="Start Analysis" isSaving={isSaving} />;
  };

  setSaving = (isSaving, validationErrors = {}) => {
    this.setState({
      isSaving,
      validationErrors
    });
  };

  createJob = async (values): Promise<?Job> => {
    const { code, name, parameters: formParams } = values;
    const {
      stats,
      filters: { enabled }
    } = formParams;
    const { pushNotification } = this.props;
    const parameters = {
      source_sample_group: values.source_sample_group,
      sample_type: values.sample_type,
      condition_variables: values.condition_variables,
      contrasts: values.contrasts,
      parameters: {
        pcut: formParams.pcut,
        log_offset: formParams.log_offset,
        when_apply_filter: formParams.when_apply_filter,
        norm: formParams.norm,
        norm_args: {
          method: formParams.norm_args.method,
          locfunc: formParams.norm_args.locfunc
        },
        stats,
        stats_args: {
          deseq: stats.includes('deseq') ? formParams.stats_args.deseq : null,
          edger: stats.includes('edger') ? formParams.stats_args.edger : null,
          limma: stats.includes('limma') ? formParams.stats_args.limma : null
        },
        filters: {
          length: enabled.includes('length')
            ? {
                length: formParams.filters.length.length
              }
            : null,
          avg_reads: enabled.includes('reads')
            ? {
                average_per_bp: formParams.filters.avg_reads.average_per_bp,
                quantile: formParams.filters.avg_reads.quantile
              }
            : null,
          expression: enabled.includes('expression')
            ? {
                median: formParams.filters.expression.median,
                mean: formParams.filters.expression.mean,
                quantile:
                  formParams.filters.expression.quantile === ''
                    ? null
                    : formParams.filters.expression.quantile,
                known: formParams.filters.expression.known
              }
            : null,
          presence: enabled.includes('presence')
            ? {
                frac: formParams.filters.presence.frac,
                min_count: formParams.filters.presence.min_count,
                per_condition: formParams.filters.presence.per_condition
              }
            : null
        },
        adjust_method: formParams.adjust_method,
        meta_p_method: formParams.meta_p_method,
        fig_formats: formParams.fig_formats,
        num_cores: formParams.num_cores
      }
    };
    const data = await Api.Analysis.createDiffExpAnalysis(
      code,
      name,
      parameters
    );
    if (data.validationErrors) {
      pushNotification(
        'Errors occurred during validation of input parameters. Please review the form!',
        'warning'
      );
      this.setSaving(false, data.validationErrors);
      return null;
    }
    const { data: job } = data;
    pushNotification('Sample group created!');
    return job;
  };

  formSubmit = async values => {
    const { pushNotification, redirect, refreshJobs } = this.props;
    this.setSaving(true);
    try {
      console.log(values);
      const job = await this.createJob(values);
      if (job) {
        await Api.Jobs.submitJob(job.id);
        refreshJobs();
        redirect(JOBS);
      }
    } catch (e) {
      pushNotification(`An error occurred: ${e.message}`, 'error');
      this.setSaving(false);
    }
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
                    <div>{this.getStep3(values)}</div>
                    <div>{this.getStep4(values)}</div>
                    <div>{this.getStep5(values)}</div>
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
