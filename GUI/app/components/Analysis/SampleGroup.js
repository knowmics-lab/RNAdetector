// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import Paper from '@material-ui/core/Paper';
import Box from '@material-ui/core/Box';
import FormGroup from '@material-ui/core/FormGroup';
import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import CircularProgress from '@material-ui/core/CircularProgress';
import { green } from '@material-ui/core/colors';
import { Formik, Form } from 'formik';
import * as Yup from 'yup';
import Backdrop from '@material-ui/core/Backdrop';
import { InputLabel } from '@material-ui/core';
import * as Api from '../../api';
import { JOBS } from '../../constants/routes';
import SelectField from '../Form/SelectField';
import TextField from '../Form/TextField';
import Wizard from '../UI/Wizard';
import FileSelector from '../UI/FileSelector';
import type { File } from '../UI/FileSelector';
import UploadProgress from '../UI/UploadProgress';
import SwitchField from '../Form/SwitchField';
import type { Job } from '../../types/jobs';
import TableField from '../Form/TableField';

type Props = {
  refreshJobs: () => void,
  redirect: mixed => void,
  pushNotification: (
    string,
    ?('success' | 'warning' | 'error' | 'info')
  ) => void,
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
  buttonWrapper: {
    margin: theme.spacing(1),
    position: 'relative'
  },
  buttonProgress: {
    color: green[500],
    position: 'absolute',
    top: '50%',
    left: '50%',
    marginTop: -12,
    marginLeft: -12
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
  descriptionFile: ?File,
  validationErrors: *,
  isUploading: boolean,
  uploadFile: string,
  uploadedBytes: number,
  uploadedPercent: number,
  uploadTotal: number
};

class SampleGroup extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.state = {
      isLoading: false,
      isSaving: false,
      ...Api.Upload.ui.initUploadState(),
      descriptionFile: null,
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
      de_novo: Yup.boolean().notRequired(),
      analysis: Yup.string()
        .required()
        .oneOf(Object.keys(SUPPORTED_ANALYSIS)),
      jobs: Yup.array()
        .required()
        .min(1)
        .of(Yup.number())
    });

  getSteps = () => ['Choose name and type', 'Select analysis', 'Upload files'];

  getStep0 = () => {
    const { classes } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can choose group code, group name, and type of contained
          samples. If you with to create a new annotation for the samples, you
          can select the &quot;De Novo Annotation&quot; option. This option will
          create a novel samples annotation using the provided description file.
          Sample code is used for further analysis, such as Differential
          Expression Analysis. It identifies the sample during the analysis
          therefore it should be a string without any spaces (only letters,
          numbers, and dashes).
        </Typography>
        <TextField label="Sample Code" name="code" required />
        <TextField label="Analysis Name" name="name" required />
        <SelectField
          label="Samples Type"
          name="analysis"
          options={SUPPORTED_ANALYSIS}
          required
        />
        <SwitchField label="De Novo Annotation?" name="de_novo" />
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

  handleDescriptionFileAdd = (descriptionFile: File[]) =>
    this.setState({
      descriptionFile: descriptionFile[0]
    });

  handleDescriptionFileRemove = () =>
    this.setState({
      descriptionFile: null
    });

  getStep2 = () => {
    const { classes } = this.props;
    const {
      isUploading,
      uploadFile,
      uploadedBytes,
      uploadedPercent,
      uploadTotal
    } = this.state;
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
              <FileSelector
                onFileRemove={this.handleDescriptionFileRemove}
                onFileAdd={this.handleDescriptionFileAdd}
                filters={[
                  {
                    name: 'TSV files',
                    extensions: ['tab', 'tsv', 'txt']
                  }
                ]}
                disabled={isUploading}
              />
            </Grid>
          </Grid>
        </FormGroup>
        <UploadProgress
          isUploading={isUploading}
          uploadFile={uploadFile}
          uploadedBytes={uploadedBytes}
          uploadedPercent={uploadedPercent}
          uploadTotal={uploadTotal}
        />
      </>
    );
  };

  getSubmitButton = () => {
    const { classes } = this.props;
    const { isSaving } = this.state;
    return (
      <>
        <Button
          type="submit"
          variant="contained"
          color="primary"
          disabled={isSaving}
        >
          Save
        </Button>
        {isSaving && (
          <CircularProgress size={24} className={classes.buttonProgress} />
        )}
      </>
    );
  };

  setSaving = (isSaving, validationErrors = {}) => {
    this.setState({
      isSaving,
      validationErrors
    });
  };

  uploadFile = async (job: Job, file: File) => {
    Api.Upload.ui.uploadStart(this.setState.bind(this), file.name);
    await Api.Upload.upload(
      job,
      file.path,
      file.name,
      file.type,
      Api.Upload.ui.makeOnProgress(this.setState.bind(this))
    );
    Api.Upload.ui.uploadEnd(this.setState.bind(this));
  };

  createGroup = async (values): Promise<?Job> => {
    // eslint-disable-next-line camelcase
    const { code, name, de_novo, jobs } = values;
    const { descriptionFile } = this.state;
    const { pushNotification } = this.props;
    const data = await Api.Analysis.createSampleGroup(
      code,
      name,
      jobs,
      descriptionFile ? descriptionFile.name : undefined,
      de_novo
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
    if (descriptionFile) {
      await this.uploadFile(job, descriptionFile);
    }
    pushNotification('Sample group created!');
    return job;
  };

  formSubmit = async values => {
    const { pushNotification, redirect, refreshJobs } = this.props;
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
              Long RNAs Analysis
            </Typography>
            <Formik
              initialValues={{
                code: '',
                name: '',
                de_novo: false,
                analysis: Object.keys(SUPPORTED_ANALYSIS)[1],
                jobs: []
              }}
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

export default withStyles(style)(SampleGroup);
