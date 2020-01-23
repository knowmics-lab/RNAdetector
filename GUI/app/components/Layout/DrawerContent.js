/* eslint-disable react/jsx-fragments */
// @flow
import * as React from 'react';
import { withStyles } from '@material-ui/core/styles';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import ListItemText from '@material-ui/core/ListItemText';
import Icon from '@material-ui/core/Icon';
import { NavLink as RouterLink } from 'react-router-dom';
import ExpandLess from '@material-ui/icons/ExpandLess';
import ExpandMore from '@material-ui/icons/ExpandMore';
import Collapse from '@material-ui/core/Collapse';
import { items as menuItems } from '../../constants/menu.json';
import * as Api from '../../api';

type ListItemLinkProps = {
  icon?: ?React.Element<*>,
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
  icon?: ?React.Element<*>,
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

export type MenuItem = {
  icon: string,
  text: string,
  collapsible: boolean,
  configured: boolean,
  key: string,
  items?: MenuItem[],
  to?: string
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
    return collapsibleState[panelKey] || false;
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

  renderMenuItems = (
    items: ?(MenuItem[]),
    nested: boolean = false
  ): ?React.Element<*> => {
    if (!items) return null;
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
    { icon, text, collapsible, configured, key, items, to }: MenuItem,
    nested: boolean = false
  ): ?React.Element<*> => {
    if (configured && !Api.Settings.isConfigured()) return null;
    const { classes } = this.props;
    if (collapsible) {
      return (
        <React.Fragment key={key}>
          <ListItemExpandable
            icon={<Icon className={icon} />}
            primary={text}
            isOpen={this.getCollapsibleState(key)}
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
    if (to) {
      return (
        <ListItemLink
          icon={<Icon className={icon} />}
          primary={text}
          to={to}
          key={key}
          className={nested ? classes.nested : null}
        />
      );
    }
    return null;
  };

  render() {
    return <>{this.renderMenuItems(menuItems)}</>;
  }
}

export default withStyles(style)(DrawerContent);
