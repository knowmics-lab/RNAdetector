/* eslint-disable react/no-unused-prop-types,camelcase */
// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import Paper from '@material-ui/core/Paper';
import Box from '@material-ui/core/Box';
import Grid from '@material-ui/core/Grid';
import CircularProgress from '@material-ui/core/CircularProgress';
import { Formik, Form } from 'formik';
import * as Yup from 'yup';
import Backdrop from '@material-ui/core/Backdrop';
import { api } from 'electron-util';
import { Collapse } from '@material-ui/core';
import { DockerManager, resetInstance } from '../../api/docker';
import * as Api from '../../api';
import SelectField from '../Form/SelectField';
import TextField from '../Form/TextField';
import Wizard from '../UI/Wizard';
import SwitchField from '../Form/SwitchField';
import { SubmitButton } from '../UI/Button';
import FileField from '../Form/FileField';
import Check from '../../api/check';
import type { ConfigObjectType } from '../../types/settings';

type Props = {
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
  freePort: number,
  logContent: string
};

class SetupWizard extends React.Component<Props, State> {
  props: Props;

  cachedContrasts: ?{
    key: ?string,
    contrasts: ?{ [string]: string }
  };

  static settingsFromValues(values): ConfigObjectType {
    return {
      local: values.local,
      apiProtocol: values.apiProtocol,
      apiHostname: values.apiHostname,
      apiPort: values.apiPort,
      apiPath: values.apiPath,
      publicPath: values.publicPath,
      dataPath: values.dataPath,
      socketPath: values.socketPath,
      containerName: values.containerName,
      apiKey: values.apiKey
    };
  }

  constructor(props) {
    super(props);
    this.cachedContrasts = null;
    this.state = {
      isLoading: true,
      isSaving: false,
      freePort: 9898,
      logContent: ''
    };
  }

  componentDidMount(): void {
    Api.Settings.findFreePort(9898)
      .then(freePort =>
        this.setState({
          isLoading: false,
          freePort
        })
      )
      .catch(e => console.error(e));
  }

  // noinspection JSCheckFunctionSignatures
  componentDidUpdate() {
    if (this.logRef.current) {
      this.logRef.current.scrollIntoView({ behavior: 'smooth' });
    }
  }

  logRef = React.createRef();

  getValidationSchema = () =>
    Yup.object().shape({
      local: Yup.boolean(),
      dataPath: Yup.string().when('local', {
        is: true,
        then: Yup.string().required(),
        otherwise: Yup.string().notRequired()
      }),
      socketPath: Yup.string().notRequired(),
      containerName: Yup.string().when('local', {
        is: true,
        then: Yup.string().required(),
        otherwise: Yup.string().notRequired()
      }),
      apiProtocol: Yup.string().required(),
      apiHostname: Yup.string().required(),
      apiPort: Yup.number()
        .positive()
        .integer()
        .min(1)
        .max(65535),
      apiPath: Yup.string().required(),
      publicPath: Yup.string().required(),
      apiKey: Yup.string().when('local', {
        is: true,
        then: Yup.string().notRequired(),
        otherwise: Yup.string().required()
      })
    });

  getSteps = () => [
    'Docker parameters',
    'Connection parameters',
    'Complete setup'
  ];

  getStep0 = values => {
    const { classes } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can select whether you wish to use RNAdetector with a local
          docker installation or with a remote server.
        </Typography>
        <SwitchField label="Use a local docker installation?" name="local" />
        <Collapse in={values.local}>
          <FileField
            label="Local storage path"
            name="dataPath"
            dialogOptions={{ properties: ['openDirectory'] }}
            helperText="Path where all data files will be stored"
          />
          <TextField label="Local container name" name="containerName" />
          <FileField
            label="Local docker socket"
            name="socketPath"
            dialogOptions={{ properties: ['openFile'], filters: [] }}
            helperText="Leave empty for automatic detection"
          />
        </Collapse>
      </>
    );
  };

  getStep1 = values => {
    const { classes } = this.props;
    const { local } = values;
    return (
      <>
        <Typography className={classes.instructions}>
          Here you can configure the connection with RNAdetector
          {local ? ' local docker container' : ' remote server'}.
        </Typography>
        <Grid container justify="space-around" alignItems="center" spacing={3}>
          <Grid item xs>
            <SelectField
              label="API Protocol"
              name="apiProtocol"
              options={{ http: 'http', https: 'https' }}
              required
            />
          </Grid>
          <Grid item xs>
            <TextField label="API Hostname" name="apiHostname" required />
          </Grid>
          <Grid item xs>
            <TextField label="API Port" name="apiPort" type="number" required />
          </Grid>
        </Grid>
        <TextField label="API Path" name="apiPath" required />
        <TextField label="Public Path" name="publicPath" required />
        {!local && <TextField label="API key" name="apiKey" />}
      </>
    );
  };

  getStep2 = () => {
    const { classes } = this.props;
    return (
      <>
        <Typography className={classes.instructions}>
          All parameters have been set. Click &quot;Install&quot; to start the
          process.
        </Typography>
      </>
    );
  };

  getSubmitButton = () => {
    const { isSaving } = this.state;
    return <SubmitButton text="Install" isSaving={isSaving} />;
  };

  setSaving = isSaving => {
    this.setState({
      isSaving
    });
  };

  setLog = logContent =>
    this.setState({
      logContent
    });

  formSubmit = async values => {
    this.setSaving(true);
    const newSettings = SetupWizard.settingsFromValues(values);
    let log = '';
    try {
      if (newSettings.local) {
        const manager = new DockerManager(newSettings);
        if (!(await manager.hasImage())) {
          log += 'Container image not found...Downloading...\n';
          this.setLog(log);
          const state = await manager.pullImage(partialState => {
            this.setLog(`${log}${partialState.toString()}`);
          });
          log += `${state.toString()}\n`;
        }
        log += 'Validating configuration:\n';
        this.setLog(log);
        const checkedSettings = await Check.checkConfig(
          newSettings,
          Api.Settings.getConfig(),
          msg => {
            log += msg;
            this.setLog(log);
          }
        );
        log += 'Saving configuration...';
        this.setLog(log);
        // Api.Settings.saveConfig(checkedSettings);
        log += 'Ok!\n';
        this.setLog(log);
        resetInstance();
      }
    } catch (e) {
      log += `An error occurred: ${e.message}\n`;
      this.setLog(log);
    }
  };

  render() {
    const { classes } = this.props;
    const { isLoading, freePort, isSaving, logContent } = this.state;
    const steps = this.getSteps();
    return (
      <>
        <Box>
          <Paper className={classes.root}>
            {!isSaving ? (
              <>
                <Typography variant="h5" component="h3">
                  Setup Wizard
                </Typography>
                <Formik
                  initialValues={{
                    local: true,
                    apiProtocol: 'http',
                    apiHostname: 'localhost',
                    apiPort: freePort,
                    apiPath: '/api/',
                    publicPath: '/storage/',
                    dataPath: `${api.app.getPath('home')}/.RNADetector`,
                    socketPath: '',
                    containerName: 'RNAdetector',
                    apiKey: ''
                  }}
                  validationSchema={this.getValidationSchema()}
                  onSubmit={v => {
                    this.formSubmit(v).catch(() => false);
                  }}
                >
                  {({ values }) => (
                    <Form>
                      <Wizard steps={steps} submitButton={this.getSubmitButton}>
                        <div>{this.getStep0(values)}</div>
                        <div>{this.getStep1(values)}</div>
                        <div>{this.getStep2()}</div>
                      </Wizard>
                    </Form>
                  )}
                </Formik>
              </>
            ) : (
              <>
                <Typography variant="h5" component="h3">
                  Installing...
                </Typography>
                <pre>{logContent}</pre>
                <div ref={this.logRef} />
              </>
            )}
          </Paper>
        </Box>
        <Backdrop className={classes.backdrop} open={isLoading}>
          <CircularProgress color="inherit" />
        </Backdrop>
      </>
    );
  }
}

export default withStyles(style)(SetupWizard);
