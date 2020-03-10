/* eslint-disable react/no-unused-prop-types,no-unused-vars */
// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import Paper from '@material-ui/core/Paper';
import Box from '@material-ui/core/Box';
import FormGroup from '@material-ui/core/FormGroup';
import Grid from '@material-ui/core/Grid';
import CircularProgress from '@material-ui/core/CircularProgress';
import { Formik, Form } from 'formik';
import * as Yup from 'yup';
import Backdrop from '@material-ui/core/Backdrop';
import { InputLabel } from '@material-ui/core';
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
  selectedJob: ?Job
};

class DiffExpr extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.state = {
      isLoading: false,
      isSaving: false,
      validationErrors: {},
      selectedJob: null
    };
  }

  loadJob = id => {
    if (!id) {
      this.setState({
        selectedJob: null
      });
    } else {
      this.setState({
        isLoading: true
      });
      const { pushNotification } = this.props;
      Api.Jobs.fetchJobById(id)
        .then(selectedJob =>
          this.setState({
            isLoading: false,
            selectedJob
          })
        )
        .catch(e => {
          pushNotification(`An error occurred: ${e.message}!`, 'error');
          this.setState({
            isLoading: false
          });
        });
    }
  };

  getVariables() {
    const { selectedJob } = this.state;
    if (
      !selectedJob ||
      !selectedJob.output ||
      !Array.isArray(selectedJob.output.metadata)
    ) {
      return [];
    }
    return selectedJob.output.metadata
      .map(f => (f.type === 'string' ? f.field : null))
      .filter(v => v !== null);
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
          name="jobs"
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

  getStep1 = values => {
    const { classes, pushNotification } = this.props;
    const { analysis } = values;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can choose which samples will be included in this group. If
          you select another sample group, all its samples and annotations will
          be included (Annotations are ignored if De Novo is selected).
        </Typography>
        <TableField
          name="jobs"
          required
          getData={() => Api.Jobs.fetchAllByType(analysis)}
          onError={e =>
            pushNotification(`An error occurred: ${e.message}!`, 'error')
          }
          label="Select one or more analysis"
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
              dataField: 'readable_type',
              sortingField: 'job_type',
              label: 'Type'
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
        <FormGroup row className={classes.formControl}>
          <Grid container alignItems="center" spacing={3}>
            <Grid item xs="auto">
              <InputLabel>Select the optional description table</InputLabel>
            </Grid>
            <Grid item xs>
              TODO
            </Grid>
          </Grid>
        </FormGroup>
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
              initialValues={DiffExpDefaults}
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
