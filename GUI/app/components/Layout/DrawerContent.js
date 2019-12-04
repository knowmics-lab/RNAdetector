// @flow
import * as React from 'react';
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

type ListItemLinkProps = {
  icon?: ?JSX.Element,
  primary: string,
  to: string
};

const ListItemLink = ({ icon, primary, to }: ListItemLinkProps) => {
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
    <li>
      <ListItem button component={renderLink}>
        {icon ? <ListItemIcon>{icon}</ListItemIcon> : null}
        <ListItemText primary={primary} />
      </ListItem>
    </li>
  );
};

ListItemLink.defaultProps = {
  icon: null
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
  console.log(isOpen);
  return (
    <li>
      <ListItem button onClick={handleClick}>
        {icon ? <ListItemIcon>{icon}</ListItemIcon> : null}
        <ListItemText primary={isOpen ? 'OPEN' : 'CLOSED'} />
        {isOpen ? <ExpandLess /> : <ExpandMore />}
      </ListItem>
    </li>
  );
};

ListItemExpandable.defaultProps = {
  icon: null
};

export default class DrawerContent extends React.Component<*, *> {
  constructor(props: *) {
    super(props);

    this.state = {
      analysisOpen: true
    };
  }

  handleAnalysisOpen = () => {
    // eslint-disable-next-line react/destructuring-assignment
    this.state.analysisOpen = !this.state.analysisOpen;
  };

  render() {
    const { analysisOpen } = this.state;
    console.log('RENDERING');
    return (
      <>
        <List>
          <ListItemExpandable
            primary="Analysis"
            isOpen={analysisOpen}
            handleClick={this.handleAnalysisOpen}
          />
          <ListItemLink
            icon={<GestureIcon />}
            primary="SmallRNA Analysis"
            to="/pippo"
          />
          <ListItemLink
            icon={<GestureIcon />}
            primary="LongRNA Analysis"
            to="/pluto"
          />
          <ListItemLink
            icon={<GestureIcon />}
            primary="CircRNA Analysis"
            to="/paperino"
          />
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
