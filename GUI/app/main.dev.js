/* eslint-disable promise/catch-or-return */
/* eslint-disable promise/always-return */
/* eslint global-require: off */

/**
 * This module executes inside of electron's main process. You can start
 * electron renderer process from here and communicate with the other processes
 * through IPC.
 *
 * When running `yarn build` or `yarn build-main`, this file is compiled to
 * `./app/main.prod.js` using webpack. This gives us some performance wins.
 *
 * @flow
 */
import { app, BrowserWindow, ipcMain, dialog } from 'electron';
import { autoUpdater } from 'electron-updater';
import log from 'electron-log';
import { download } from 'electron-dl';
import tus from 'tus-js-client';
import fs from 'fs';
import MenuBuilder from './menu';
import { Settings, Docker } from './api';

export default class AppUpdater {
  constructor() {
    log.transports.file.level = 'info';
    autoUpdater.logger = log;
    autoUpdater.checkForUpdatesAndNotify();
  }
}

let mainWindow = null;

if (process.env.NODE_ENV === 'production') {
  const sourceMapSupport = require('source-map-support');
  sourceMapSupport.install();
}

if (
  process.env.NODE_ENV === 'development' ||
  process.env.DEBUG_PROD === 'true'
) {
  require('electron-debug')();
}

const installExtensions = async () => {
  const installer = require('electron-devtools-installer');
  const forceDownload = !!process.env.UPGRADE_EXTENSIONS;
  const extensions = ['REACT_DEVELOPER_TOOLS', 'REDUX_DEVTOOLS'];

  return Promise.all(
    extensions.map(name => installer.default(installer[name], forceDownload))
  ).catch(console.log);
};

/**
 * Add event listeners...
 */

app.on('window-all-closed', () => {
  // Respect the OSX convention of having the application in memory even
  // after all windows have been closed
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('ready', async () => {
  if (
    process.env.NODE_ENV === 'development' ||
    process.env.DEBUG_PROD === 'true'
  ) {
    await installExtensions();
  }

  mainWindow = new BrowserWindow({
    show: false,
    width: 1024,
    height: 728,
    webPreferences: {
      nodeIntegration: true,
      nativeWindowOpen: true
    }
  });

  mainWindow.loadURL(`file://${__dirname}/app.html`);

  // @TODO: Use 'ready-to-show' event
  //        https://github.com/electron/electron/blob/master/docs/api/browser-window.md#using-ready-to-show-event
  mainWindow.webContents.on('did-finish-load', () => {
    if (!mainWindow) {
      throw new Error('"mainWindow" is not defined');
    }
    if (process.env.START_MINIMIZED) {
      mainWindow.minimize();
    } else {
      mainWindow.show();
      mainWindow.focus();
    }
  });

  let doQuit = false;

  const quitNow = () => {
    doQuit = true;
    if (mainWindow) mainWindow.close();
    app.quit();
  };

  const stopDockerAndQuit = () => {
    console.log('Waiting for docker container to stop');
    Docker.stopContainer()
      .then(() => {
        console.log('Docker container has been stopped! Quitting!');
      })
      .catch(ex => {
        console.log('Docker container cannot be stopped! Stop it manually!');
        console.log(ex);
      })
      .finally(() => {
        quitNow();
      });
  };

  mainWindow.on('close', e => {
    if (Settings.isLocal() && Settings.isConfigured() && !doQuit) {
      e.preventDefault();
      if (Settings.autoStopDockerOnClose()) {
        stopDockerAndQuit();
      } else {
        dialog.showMessageBox(
          mainWindow,
          {
            // title: 'Close docker',
            message: 'Do you wish to close docker?',
            buttons: ['&Yes', '&No'],
            type: 'question',
            checkboxLabel: 'Close without asking again?',
            checkboxChecked: Settings.autoStopDockerOnClose(),
            normalizeAccessKeys: true
          },
          (response, checkboxChecked) => {
            if (checkboxChecked) {
              Settings.setAutoStopDockerOnClose();
            }
            if (response === 0) {
              if (mainWindow)
                mainWindow.webContents.send('waiting-docker-close');
              stopDockerAndQuit();
            } else {
              doQuit = true;
              if (mainWindow) mainWindow.close();
              quitNow();
            }
          }
        );
      }
    }
  });

  mainWindow.on('closed', async () => {
    mainWindow = null;
  });

  const menuBuilder = new MenuBuilder(mainWindow);
  menuBuilder.buildMenu();

  ipcMain.on('download-file', async (event, { id, url, filename }) => {
    const win = BrowserWindow.getFocusedWindow();
    await download(win, url, {
      saveAs: true,
      openFolderWhenDone: false,
      filename,
      onStarted() {
        event.reply('download-started', { id });
      },
      onProgress(progress) {
        if (progress.percent >= 1) {
          event.reply('download-completed', { id });
        }
      }
    });
  });

  ipcMain.on(
    'upload-file',
    async (event, { id, filePath, fileName, fileType, endpoint }) => {
      const file = fs.createReadStream(filePath);
      const { size } = fs.statSync(filePath);
      const upload = new tus.Upload(file, {
        endpoint,
        retryDelays: [0, 3000, 5000, 10000, 20000],
        headers: {
          ...Settings.getAuthHeaders()
        },
        chunkSize: 50 * 1024 * 1024, // 50Mb per chunk
        resume: true,
        metadata: {
          filename: fileName,
          filetype: fileType
        },
        uploadSize: size,
        onError(error) {
          event.reply('upload-message', {
            id,
            isDone: false,
            isProgress: false,
            isError: true,
            percentage: 0,
            bytesUploaded: 0,
            bytesTotal: 0,
            error: error.message
          });
        },
        onProgress(bytesUploaded, bytesTotal) {
          const percentage = Math.round((bytesUploaded / bytesTotal) * 100);
          event.reply('upload-message', {
            id,
            isDone: false,
            isProgress: true,
            isError: false,
            error: null,
            percentage,
            bytesUploaded,
            bytesTotal
          });
        },
        onSuccess() {
          event.reply('upload-message', {
            id,
            isDone: true,
            isProgress: false,
            isError: false,
            error: null,
            percentage: 0,
            bytesUploaded: 0,
            bytesTotal: 0
          });
        }
      });
      upload.start();
    }
  );

  mainWindow.webContents.on(
    'new-window',
    (event, url, frameName, disposition, options) => {
      event.preventDefault();
      const forcedOptions = {
        parent: mainWindow,
        width: 1024,
        height: 728,
        webPreferences: {
          nodeIntegration: false,
          nativeWindowOpen: true,
          webviewTag: false,
          nodeIntegrationInSubFrames: false
        }
      };
      const isModal = frameName === 'modal';
      const win = new BrowserWindow({
        ...options,
        ...forcedOptions,
        modal: isModal
      });
      win.setMenuBarVisibility(false);
      // eslint-disable-next-line no-param-reassign
      event.newGuest = win;
    }
  );

  // Remove this if your app does not use auto updates
  // eslint-disable-next-line
  new AppUpdater();
});
