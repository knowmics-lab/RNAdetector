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
import Breadcrumbs from '@material-ui/core/Breadcrumbs';
import Link from '@material-ui/core/Link';
import * as Api from '../api';
import { JOBS, REFERENCES } from '../constants/routes.json';
import SelectField from './Form/SelectField';
import TextField from './Form/TextField';
import Wizard from './UI/Wizard';
import FileSelector from './UI/FileSelector';
import type { File } from './UI/FileSelector';
import UploadProgress from './UI/UploadProgress';
import CustomArgumentsField, {
  ProcessCustomArguments
} from './Form/CustomArgumentsField';

type Props = {
  refreshJobs: () => void,
  refreshReferences: () => void,
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
    instructions: *
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
  }
});

type State = {
  isSaving: boolean,
  isUploading: boolean,
  uploadFile: string,
  uploadedBytes: number,
  uploadedPercent: number,
  uploadTotal: number,
  fastaFile: ?File,
  mapFile: ?File,
  validationErrors: *
};

class CreateReference extends React.Component<Props, State> {
  props: Props;

  constructor(props) {
    super(props);
    this.state = {
      isSaving: false,
      ...Api.Upload.ui.initUploadState(),
      fastaFile: null,
      mapFile: null,
      validationErrors: {}
    };
  }

  getValidationSchema = () =>
    Yup.object().shape({
      name: Yup.string()
        .required()
        .matches(/^[A-Za-z0-9\-_]+$/, {
          message: 'The field must contain only letters, numbers, and dashes.'
        }),
      availableFor: Yup.array()
        .of(Yup.string().oneOf(['bwa', 'hisat', 'salmon', 'star']))
        .required()
    });

  getSteps = () => ['Choose a name', 'Select aligners', 'Select a file'];

  getStep0 = () => {
    const { classes } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Choose a name fo the new reference sequence. The name should contain
          only letters, numbers, and dashes.
        </Typography>
        <TextField label="Name" name="name" required />
      </>
    );
  };

  getStep1 = values => {
    const { classes } = this.props;
    const { availableFor } = values;
    return (
      <>
        <Typography className={classes.instructions}>
          Choose which alignment algorithms will be available for this sequence.
          Based on the selected algorithms, the appropriate indexing methods
          will be used.
        </Typography>
        <SelectField
          label="Indexing algorithms"
          name="availableFor"
          options={{
            bwa: 'BWA',
            hisat: 'HISAT2',
            salmon: 'Salmon',
            star: 'STAR'
          }}
          multiple
          required
        />
        {availableFor.includes('bwa') && (
          <CustomArgumentsField
            name="custom_arguments.bwa"
            labelEnable="Enable Custom Arguments for BWA indexing?"
          />
        )}
        {availableFor.includes('hisat') && (
          <CustomArgumentsField
            name="custom_arguments.hisat"
            labelEnable="Enable Custom Arguments for HISAT2 indexing?"
          />
        )}
        {availableFor.includes('salmon') && (
          <CustomArgumentsField
            name="custom_arguments.salmon"
            labelEnable="Enable Custom Arguments for Salmon indexing?"
          />
        )}
        {availableFor.includes('star') && (
          <CustomArgumentsField
            name="custom_arguments.star"
            labelEnable="Enable Custom Arguments for STAR indexing?"
          />
        )}
      </>
    );
  };

  handleFileAdd = (field: string) => f =>
    this.setState({
      [field]: f[0]
    });

  handleFileRemove = (field: string) => () =>
    this.setState({
      [field]: null
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
          Select the FASTA file you wish to upload and click &quot;Save&quot; to
          start the upload process. If you are uploading a transcriptome and you
          wish to enable pathway analysis support, please upload also a TSV file
          mapping transcripts and gene ids to Entrez Ids.
        </Typography>
        <FormGroup row className={classes.formControl}>
          <Grid container justify="center" alignItems="flex-start">
            <Grid item xs>
              <FileSelector
                title="Select FASTA file"
                onFileRemove={this.handleFileRemove('fastaFile')}
                onFileAdd={this.handleFileAdd('fastaFile')}
                filters={[
                  { name: 'FASTA files', extensions: ['fa', 'fasta', 'gz'] }
                ]}
                disabled={isUploading}
              />
            </Grid>
            <Grid item xs>
              <FileSelector
                title="Select map file"
                onFileRemove={this.handleFileRemove('mapFile')}
                onFileAdd={this.handleFileAdd('mapFile')}
                filters={[{ name: 'TSV files', extensions: ['tsv', 'txt'] }]}
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

  formSubmit = async values => {
    const {
      pushNotification,
      redirect,
      refreshJobs,
      refreshReferences
    } = this.props;
    const { availableFor, custom_arguments: customArguments } = values;
    let finalCustomArguments = {};
    // eslint-disable-next-line no-restricted-syntax
    for (const algo of availableFor) {
      finalCustomArguments = ProcessCustomArguments(
        finalCustomArguments,
        customArguments[algo],
        algo
      );
    }
    const { fastaFile, mapFile } = this.state;
    if (!fastaFile) {
      pushNotification('You must select a FASTA file.', 'error');
    } else {
      this.setSaving(true);
      try {
        const data = await Api.References.create(
          values.name,
          fastaFile.name,
          availableFor,
          mapFile ? mapFile.name : null,
          finalCustomArguments
        );
        if (data.validationErrors) {
          pushNotification(
            'Errors occurred during validation of input parameters. Please review the form!',
            'warning'
          );
          this.setSaving(false, data.validationErrors);
        } else {
          const { data: job } = data;
          pushNotification(
            'A new indexing job has been created! Uploading FASTA file...'
          );
          const files = [fastaFile];
          if (mapFile) files.push(mapFile);
          // eslint-disable-next-line no-restricted-syntax
          for (const f of files) {
            Api.Upload.ui.uploadStart(this.setState.bind(this), f.name);
            // eslint-disable-next-line no-await-in-loop
            await Api.Upload.upload(
              job,
              f.path,
              f.name,
              f.type,
              Api.Upload.ui.makeOnProgress(this.setState.bind(this))
            );
          }
          Api.Upload.ui.uploadEnd(this.setState.bind(this));
          pushNotification('FASTA file uploaded! Starting indexing job...');
          await Api.Jobs.submitJob(job.id);
          pushNotification('Indexing job queued!');
          refreshJobs();
          refreshReferences();
          redirect(JOBS);
        }
      } catch (e) {
        pushNotification(`An error occurred: ${e.message}`, 'error');
        this.setSaving(false);
      }
    }
  };

  redirectReference = event => {
    const { redirect } = this.props;
    event.preventDefault();
    redirect(REFERENCES);
  };

  render() {
    const { classes } = this.props;
    const { validationErrors } = this.state;
    const steps = this.getSteps();
    return (
      <>
        <Breadcrumbs aria-label="breadcrumb">
          <Link color="inherit" href="#" onClick={this.redirectReference}>
            References
          </Link>
          <Typography color="textPrimary">Add Reference</Typography>
        </Breadcrumbs>
        <Box>
          <Paper className={classes.root}>
            <Typography variant="h5" component="h3">
              Add reference genome/transcriptome
            </Typography>
            <Formik
              initialValues={{
                name: '',
                availableFor: [],
                custom_arguments: {
                  bwa: { enable: false, value: '' },
                  salmon: { enable: false, value: '' },
                  hisat: { enable: false, value: '' },
                  star: { enable: false, value: '' }
                }
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
      </>
    );
  }
}

export default withStyles(style)(CreateReference);
