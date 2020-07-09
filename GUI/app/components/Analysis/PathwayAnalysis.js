/* eslint-disable react/no-unused-prop-types,camelcase */
// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import Paper from '@material-ui/core/Paper';
import Box from '@material-ui/core/Box';
import CircularProgress from '@material-ui/core/CircularProgress';
import { Formik, Form } from 'formik';
import * as Yup from 'yup';
import Backdrop from '@material-ui/core/Backdrop';
import * as Api from '../../api';
import { JOBS } from '../../constants/routes.json';
import SelectField from '../Form/SelectField';
import TextField from '../Form/TextField';
import Wizard from '../UI/Wizard';
import SwitchField from '../Form/SwitchField';
import type { Job } from '../../types/jobs';
import TableField from '../Form/TableField';
import type { PushNotificationFunction } from '../../types/notifications';
import { SubmitButton } from '../UI/Button';
import FormGroup from "@material-ui/core/FormGroup";
import Grid from "@material-ui/core/Grid";
import ValidationError from "../../errors/ValidationError";

const VALID_ORGANISMS = {
  hsa: 'Human (Homo Sapiens)',
  mmu: 'Mouse (Mus Musculus)',
  rno: 'Rat (Rattus Norvegicus)'
};

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

type State = {
  isLoading: boolean,
  isSaving: boolean,
  hasValidationErrors: boolean,
  validationErrors: *
};

class PathwayAnalysis extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.state = {
      isLoading: false,
      isSaving: false,
      hasValidationErrors: false,
      validationErrors: {}
    };
  }

  getValidationSchema = () =>
    Yup.object().shape({
      code: Yup.string()
        .required()
        .matches(/^[A-Za-z0-9\-_]+$/, {
          message: 'The field must contain only letters, numbers, and dashes.'
        }),
      name: Yup.string().required(),
      degs_analysis: Yup.number().required(),
      degs: Yup.object().shape({
        p_cutoff: Yup.number()
          .positive()
          .min(0)
          .max(1),
        p_use_fdr: Yup.boolean(),
        lfc_threshold: Yup.number().min(0)
      }),
      pathways: Yup.object().shape({
        p_cutoff: Yup.number()
          .positive()
          .min(0)
          .max(1),
        p_use_fdr: Yup.boolean(),
        organism: Yup.string().oneOf(Object.keys(VALID_ORGANISMS))
      })
    });

  getSteps = () => [
    'Choose name and type',
    'DEGs Analysis Parameters',
    'Pathway Analysis Parameters'
  ];

  getStep0 = () => {
    const { classes, pushNotification } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can choose code and name for the analysis, and select a DEGs
          analysis that will be used to gather data and contrasts. Sample code
          is used for further analysis. It identifies the sample during the
          analysis therefore it should be a string without any spaces (only
          letters, numbers, and dashes).
        </Typography>
        <TextField label="Code" name="code" required />
        <TextField label="Analysis Name" name="name" required />
        <TableField
          name="degs_analysis"
          required
          single
          getData={() => Api.Jobs.fetchAllByType('diff_expr_analysis_job_type')}
          onError={e =>
            pushNotification(`An error occurred: ${e.message}!`, 'error')
          }
          label="Select one DEGs analysis"
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

  getStep1 = () => {
    const { classes } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can select the common parameters to extract DEGs.
        </Typography>
        <TextField
          label="p-value Threshold"
          name="degs.p_cutoff"
          required
          type="number"
          helperText="A p-value cutoff for exporting differentially expressed genes."
        />
        <SwitchField
          label="Should p-value cutoff be applied to FDR p-values?"
          name="degs.p_use_fdr"
          required
        />
        <TextField
          label="Log-Fold-Change Threshold"
          name="degs.lfc_threshold"
          required
          type="number"
          helperText="A minimum absolute Log-Fold-Change cutoff for exporting differentially expressed genes."
        />
      </>
    );
  };

  getStep2 = () => {
    const { classes } = this.props;
    const { hasValidationErrors, validationErrors } = this.state;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can select pathway analysis parameters. After entering all
          the parameters, you can start the analysis by clicking on the
          &quot;Start Analysis&quot; button.
        </Typography>
        <SelectField
          label="Which organism should be used?"
          name="pathways.organism"
          options={VALID_ORGANISMS}
          required
        />
        <TextField
          label="p-value Threshold"
          name="pathways.p_cutoff"
          required
          type="number"
          helperText="A p-value cutoff for exporting significantly impacted pathways."
        />
        <SwitchField
          label="Should p-value cutoff be applied to FDR p-values?"
          name="pathways.p_use_fdr"
          required
        />
        {hasValidationErrors && (
          <FormGroup row className={classes.formControl}>
            <Grid container alignItems="center" spacing={3}>
              <Grid item xs="auto">
                <Typography color="error" paragraph variant="caption">
                  Error log:
                </Typography>
                <pre>{JSON.stringify(validationErrors)}</pre>
              </Grid>
            </Grid>
          </FormGroup>
        )}
      </>
    );
  };

  getSubmitButton = () => {
    const { isSaving } = this.state;
    return <SubmitButton text="Start Analysis" isSaving={isSaving} />;
  };

  setSaving = (isSaving, validationErrors = {}) => {
    const hasValidationErrors = !(
      Object.keys(validationErrors).length === 0 &&
      validationErrors.constructor === Object
    );
    this.setState({
      isSaving,
      hasValidationErrors,
      validationErrors
    });
  };

  createJob = async (values): Promise<?Job> => {
    const { code, name } = values;
    const { pushNotification } = this.props;
    const parameters = {
      degs_analysis: values.degs_analysis,
      degs: {
        p_cutoff: values.degs.p_cutoff,
        p_use_fdr: values.degs.p_use_fdr,
        lfc_threshold: values.degs.lfc_threshold
      },
      pathways: {
        organism: values.pathways.organism,
        p_cutoff: values.pathways.p_cutoff,
        p_use_fdr: values.pathways.p_use_fdr
      }
    };
    const data = await Api.Analysis.createPathwayAnalysis(
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
    pushNotification('Pathway analysis created!');
    return job;
  };

  formSubmit = async values => {
    const { pushNotification, redirect, refreshJobs } = this.props;
    this.setSaving(true);
    try {
      const job = await this.createJob(values);
      if (job) {
        await Api.Jobs.submitJob(job.id);
        refreshJobs();
        redirect(JOBS);
      }
    } catch (e) {
      pushNotification(`An error occurred: ${e.message}`, 'error');
      if (!(e instanceof ValidationError)) {
        this.setSaving(false);
      }
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
              Pathway Analysis
            </Typography>
            <Formik
              initialValues={{
                code: '',
                name: '',
                degs_analysis: null,
                degs: {
                  p_cutoff: 0.05,
                  p_use_fdr: true,
                  lfc_threshold: 0.0
                },
                pathways: {
                  organism: 'hsa',
                  p_cutoff: 0.05,
                  p_use_fdr: true
                }
              }}
              initialErrors={validationErrors}
              validationSchema={this.getValidationSchema()}
              onSubmit={v => {
                this.formSubmit(v).catch(() => false);
              }}
            >
              <Form>
                <Wizard steps={steps} submitButton={this.getSubmitButton}>
                  <div>{this.getStep0()}</div>
                  <div>{this.getStep1()}</div>
                  <div>{this.getStep2()}</div>
                </Wizard>
              </Form>
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

export default withStyles(style)(PathwayAnalysis);
