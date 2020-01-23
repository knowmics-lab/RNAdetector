// @flow
import * as React from 'react';
import { makeStyles, Theme } from '@material-ui/core/styles';
import Drawer from '@material-ui/core/Drawer';
import AppBar from '@material-ui/core/AppBar';
import CssBaseline from '@material-ui/core/CssBaseline';
import Toolbar from '@material-ui/core/Toolbar';
import Typography from '@material-ui/core/Typography';
import Divider from '@material-ui/core/Divider';
import MenuIcon from '@material-ui/icons/Menu';
import IconButton from '@material-ui/core/IconButton';
import DrawerContent from '../components/Layout/DrawerContent';
import NotificationsList from '../components/Layout/Notifications';
import CloseHandler from '../components/Layout/CloseHandler';
import StartHandler from '../components/Layout/StartHandler';

const drawerWidth = 260;
const fixDriverMinSize = 1100;

const useStyles = makeStyles((theme: Theme) => ({
  root: {
    display: 'flex'
  },
  drawer: {
    [theme.breakpoints.up(fixDriverMinSize)]: {
      width: drawerWidth,
      flexShrink: 0
    }
  },
  hideMd: {
    [theme.breakpoints.up(fixDriverMinSize)]: {
      display: 'none'
    }
  },
  showMd: {
    [theme.breakpoints.down(fixDriverMinSize)]: {
      display: 'none'
    }
  },
  appBar: {
    zIndex: theme.zIndex.drawer + 1
  },
  menuButton: {
    marginRight: theme.spacing(2),
    [theme.breakpoints.up(fixDriverMinSize)]: {
      display: 'none'
    }
  },
  toolbarResponsive: {
    [theme.breakpoints.down(fixDriverMinSize)]: {
      display: 'none'
    },
    [theme.breakpoints.up(fixDriverMinSize)]: theme.mixins.toolbar
  },
  toolbar: theme.mixins.toolbar,
  drawerPaper: {
    width: drawerWidth
  },
  content: {
    flexGrow: 1,
    padding: theme.spacing(3)
  }
}));

type Props = {
  children: React.Node | React.Node[]
};

export default function App({ children }: Props) {
  const classes = useStyles();
  const [mobileOpen, setMobileOpen] = React.useState(false);

  const handleDrawerToggle = () => {
    setMobileOpen(!mobileOpen);
  };

  const drawer = (
    <>
      <div className={classes.toolbarResponsive} />
      <Divider />
      <DrawerContent />
    </>
  );

  return (
    <div className={classes.root}>
      <CssBaseline />
      <AppBar position="fixed" className={classes.appBar}>
        <Toolbar>
          <IconButton
            color="inherit"
            aria-label="open drawer"
            edge="start"
            onClick={handleDrawerToggle}
            className={classes.menuButton}
          >
            <MenuIcon />
          </IconButton>
          <Typography variant="h6" noWrap>
            RNAdetector
          </Typography>
        </Toolbar>
      </AppBar>
      <nav className={classes.drawer}>
        <Drawer
          variant="temporary"
          anchor="left"
          open={mobileOpen}
          onClose={handleDrawerToggle}
          className={classes.hideMd}
          classes={{
            paper: classes.drawerPaper
          }}
          ModalProps={{
            keepMounted: true
          }}
        >
          {drawer}
        </Drawer>
        <Drawer
          className={classes.showMd}
          classes={{
            paper: classes.drawerPaper
          }}
          variant="permanent"
          open
        >
          {drawer}
        </Drawer>
      </nav>
      <main className={classes.content}>
        <StartHandler />
        <div className={classes.toolbar} />
        {children}
        <NotificationsList />
        <CloseHandler />
      </main>
    </div>
  );
}
