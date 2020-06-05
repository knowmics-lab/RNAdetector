import * as React from 'react';
import { ipcRenderer } from 'electron';
import * as Api from '../../api';

export default function StartHandler() {
  const [first, setFirst] = React.useState(false);

  React.useEffect(() => {
    if (!first) {
      if (Api.Settings.isConfigured() && Api.Settings.isLocal()) {
        const sendMessage = (message, error) =>
          ipcRenderer.send('display-blocking-message', {
            message,
            error
          });
        const showLog = log => ipcRenderer.send('blocking-message-log', log);
        Api.Docker.startupSequence(sendMessage, showLog)
          // eslint-disable-next-line promise/always-return
          .then(() => ipcRenderer.send('hide-blocking-message'))
          .catch(e => sendMessage(e.message, true));
      }
      setFirst(true);
    }
  });

  return <></>;
}
