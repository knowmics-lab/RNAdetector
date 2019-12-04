// @flow
import * as React from 'react';
import ReactDOM from 'react-dom';
import { withStyles } from '@material-ui/core/styles';
import Divider from '@material-ui/core/Divider';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import ListItemText from '@material-ui/core/ListItemText';
import GestureIcon from '@material-ui/icons/Gesture';
import InboxIcon from '@material-ui/icons/MoveToInbox';
import MailIcon from '@material-ui/icons/Mail';
import { NavLink as RouterLink } from 'react-router-dom';
import ExpandLess from '@material-ui/icons/ExpandLess';
import ExpandMore from '@material-ui/icons/ExpandMore';
import Collapse from '@material-ui/core/Collapse';

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
  handleClick: () => void
};

const ListItemExpandable = ({
  icon,
  primary,
  isOpen,
  handleClick
}: ListItemExpandableProps) => {
  return (
    <div>
      <ListItem button onClick={handleClick}>
        {icon ? <ListItemIcon>{icon}</ListItemIcon> : null}
        <ListItemText primary={primary} />
        {isOpen ? <ExpandLess /> : <ExpandMore />}
      </ListItem>
    </div>
  );
};

ListItemExpandable.defaultProps = {
  icon: null
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
  analysisOpen: boolean
};

class DrawerContent extends React.Component<
  DrawerContentProps,
  DrawerContentState
> {
  constructor(props: *) {
    super(props);

    this.state = {
      analysisOpen: true
    };
  }

  handleAnalysisOpen = () => {
    this.setState(prevState => ({ analysisOpen: !prevState.analysisOpen }));
  };

  render() {
    const { classes } = this.props;
    const { analysisOpen } = this.state;
    return (
      <>
        <List component="nav" className={classes.root}>
          <ListItemExpandable
            primary="Analysis"
            isOpen={analysisOpen}
            handleClick={this.handleAnalysisOpen}
          />
          <Collapse in={analysisOpen} timeout="auto" unmountOnExit>
            <List component="div" disablePadding>
              <ListItemLink
                icon={<GestureIcon />}
                primary="SmallRNA Analysis"
                to="/pippo"
                className={classes.nested}
              />
              <ListItemLink
                icon={<GestureIcon />}
                primary="LongRNA Analysis"
                to="/pluto"
                className={classes.nested}
              />
              <ListItemLink
                icon={<GestureIcon />}
                primary="CircRNA Analysis"
                to="/paperino"
                className={classes.nested}
              />
            </List>
          </Collapse>
        </List>
        <Divider />
        <List>
          {['All mail', 'Trash', 'Spam'].map((text, index) => (
            <ListItem button key={text}>
              <ListItemIcon>
                {index % 2 === 0 ? <InboxIcon /> : <MailIcon />}
              </ListItemIcon>
              <ListItemText primary={text} />
            </ListItem>
          ))}
        </List>
      </>
    );
  }
}

export default withStyles(style)(DrawerContent);
