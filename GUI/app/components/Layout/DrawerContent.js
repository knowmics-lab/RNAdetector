/* eslint-disable react/jsx-fragments */
// @flow
import * as React from 'react';
import ReactDOM from 'react-dom';
import { withStyles } from '@material-ui/core/styles';
import Divider from '@material-ui/core/Divider';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import ListItemText from '@material-ui/core/ListItemText';
import Icon from '@material-ui/core/Icon';
import WorkIcon from '@material-ui/icons/Work';
import { NavLink as RouterLink } from 'react-router-dom';
import ExpandLess from '@material-ui/icons/ExpandLess';
import ExpandMore from '@material-ui/icons/ExpandMore';
import Collapse from '@material-ui/core/Collapse';
import { items as menuItems } from '../../constants/menu.json';

type ListItemLinkProps = {
  icon?: ?JSX.Element,
  primary: string,
  to: string,
  className?: ?string
};

const ListItemLink = ({ icon, primary, to, className }: ListItemLinkProps) => {
  const renderLink = React.useMemo(
    () =>
      React.forwardRef((itemProps, ref) => (
        <RouterLink
          to={to}
          // eslint-disable-next-line react/jsx-props-no-spreading
          {...itemProps}
          innerRef={ref}
          activeClassName="Mui-selected"
        />
      )),
    [to]
  );

  return (
    <div>
      <ListItem button component={renderLink} className={className}>
        {icon ? <ListItemIcon>{icon}</ListItemIcon> : null}
        <ListItemText primary={primary} />
      </ListItem>
    </div>
  );
};

ListItemLink.defaultProps = {
  icon: null,
  className: null
};

type ListItemExpandableProps = {
  icon?: ?JSX.Element,
  primary: string,
  isOpen: boolean,
  handleClick: () => void,
  className?: ?string
};

const ListItemExpandable = ({
  icon,
  primary,
  isOpen,
  handleClick,
  className
}: ListItemExpandableProps) => {
  return (
    <div>
      <ListItem button onClick={handleClick} className={className}>
        {icon ? <ListItemIcon>{icon}</ListItemIcon> : null}
        <ListItemText primary={primary} />
        {isOpen ? <ExpandLess /> : <ExpandMore />}
      </ListItem>
    </div>
  );
};

ListItemExpandable.defaultProps = {
  icon: null,
  className: null
};

const style = theme => ({
  root: {
    width: '100%',
    backgroundColor: theme.palette.background.paper
  },
  nested: {
    paddingLeft: theme.spacing(4)
  }
});

type DrawerContentProps = {
  classes: *
};

type DrawerContentState = {
  collapsibleState: {}
};

class DrawerContent extends React.Component<
  DrawerContentProps,
  DrawerContentState
> {
  constructor(props: *) {
    super(props);

    this.state = {
      collapsibleState: {}
    };
  }

  realGetCollapsibleState = (panelKey: string, state: DrawerContentState) => {
    const { collapsibleState } = state;
    const panelState = collapsibleState[panelKey] || false;
    return panelState;
  };

  getCollapsibleState = (panelKey: string) => {
    return this.realGetCollapsibleState(panelKey, this.state);
  };

  setCollapsibleState = (panelKey: string) => {
    this.setState(prevState => ({
      collapsibleState: {
        ...prevState.collapsibleState,
        [panelKey]: !this.realGetCollapsibleState(panelKey, prevState)
      }
    }));
  };

  getCollapsibleHandler = (panelKey: string) => {
    return () => {
      this.setCollapsibleState(panelKey);
    };
  };

  renderMenuItems = (items: mixed[], nested: boolean = false) => {
    const { classes } = this.props;
    const itemsElements = items.map(item => this.renderMenuItem(item, nested));
    if (nested) {
      return (
        <List component="div" disablePadding>
          {itemsElements}
        </List>
      );
    }

    return (
      <List component="nav" className={classes.root}>
        {itemsElements}
      </List>
    );
  };

  renderMenuItem = (
    { icon, text, collapsible, key, items, to },
    nested = false
  ) => {
    const { classes } = this.props;
    if (collapsible) {
      return (
        <React.Fragment>
          <ListItemExpandable
            icon={<Icon className={icon} />}
            primary={text}
            isOpen={this.getCollapsibleState(key)}
            key={key}
            handleClick={this.getCollapsibleHandler(key)}
            className={nested ? classes.nested : null}
          />
          <Collapse
            in={this.getCollapsibleState(key)}
            timeout="auto"
            unmountOnExit
            className={nested ? classes.nested : null}
          >
            {this.renderMenuItems(items, true)}
          </Collapse>
        </React.Fragment>
      );
    }
    return (
      <ListItemLink
        icon={<Icon className={icon} />}
        primary={text}
        to={to}
        key={key}
        className={nested ? classes.nested : null}
      />
    );
  };

  render() {
    const { classes } = this.props;
    const analysisOpen = this.getCollapsibleState('test');
    return this.renderMenuItems(menuItems);
    return (
      <>
        <List component="nav" className={classes.root}>
          <ListItemExpandable
            icon={<Icon className="far fa-play-circle" />}
            primary="Run Analysis"
            isOpen={analysisOpen}
            handleClick={this.getCollapsibleHandler('test')}
          />
          <Collapse in={analysisOpen} timeout="auto" unmountOnExit>
            <List component="div" disablePadding>
              <ListItemLink
                icon={<Icon className="fas fa-dna" />}
                primary="SmallRNA Analysis"
                to="/pippo"
                className={classes.nested}
              />
              <ListItemLink
                icon={<Icon className="fas fa-dna" />}
                primary="LongRNA Analysis"
                to="/pluto"
                className={classes.nested}
              />
              <ListItemLink
                icon={<Icon className="far fa-circle" />}
                primary="CircRNA Analysis"
                to="/paperino"
                className={classes.nested}
              />
            </List>
          </Collapse>
          <ListItemLink icon={<WorkIcon />} primary="Jobs" to="/jobs" />
        </List>
      </>
    );
  }
}

export default withStyles(style)(DrawerContent);
