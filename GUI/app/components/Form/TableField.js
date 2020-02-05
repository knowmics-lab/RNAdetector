// @flow
import * as React from 'react';
import Table from '@material-ui/core/Table';
import TablePagination from '@material-ui/core/TablePagination';
import TableContainer from '@material-ui/core/TableContainer';
import LinearProgress from '@material-ui/core/LinearProgress';
import { withStyles } from '@material-ui/core/styles';
import { connect, getIn } from 'formik';
// $FlowFixMe
import type { FormikContextType } from 'formik';
import FormControl from '@material-ui/core/FormControl';
import FormLabel from '@material-ui/core/FormLabel';
import FormHelperText from '@material-ui/core/FormHelperText';
import type { TableColumn } from '../UI/Table/types';
import TableHeader from '../UI/Table/Header';
import TableBody from '../UI/Table/Body';

export type TableProps = {
  label: string,
  required?: boolean,
  helperText?: string,
  name: string,
  size?: 'small' | 'medium',
  columns: TableColumn[],
  keyField?: string,
  getData: () => Promise<{ +[string]: * }[]>,
  formik: FormikContextType,
  onError?: (*) => void,
  classes: {
    root: *,
    container: *,
    label: *,
    stickyStyle: *,
    loading: *
  }
};

export type TableState = {
  data: ?({ +[string]: * }[]),
  currentPage: ?number,
  perPage: number,
  total: ?number,
  fetching: boolean
};

const styles = theme => ({
  root: {
    width: '100%',
    margin: theme.spacing(1),
    minWidth: 120
  },
  container: {
    // maxHeight: 440
  },
  label: {
    fontSize: '0.75rem'
  },
  stickyStyle: {
    backgroundColor: 'white'
  },
  loading: {
    width: '100%',
    '& > * + *': {
      marginTop: theme.spacing(2)
    }
  }
});

class TableField extends React.Component<TableProps, TableState> {
  props: TableProps;

  static defaultProps = {
    size: 'small',
    keyField: 'id',
    required: false,
    helperText: '',
    onError: undefined
  };

  constructor(props) {
    super(props);
    this.state = {
      data: null,
      currentPage: null,
      perPage: 10,
      total: null,
      fetching: false
    };
  }

  componentDidMount(): void {
    const { data } = this.state;
    if (data === null) {
      const { getData, onError } = this.props;
      this.setState({
        fetching: true
      });
      // eslint-disable-next-line promise/catch-or-return
      getData()
        // eslint-disable-next-line promise/always-return
        .then(d => {
          this.setState({
            data: d,
            currentPage: 0,
            total: d.length
          });
        })
        .catch(e => (typeof onError === 'function' ? onError(e) : undefined))
        .finally(() => {
          this.setState({
            fetching: false
          });
        });
    }
  }

  handleChangePage = (event, newPage) => {
    this.setState({
      currentPage: newPage
    });
  };

  handleChangeRowsPerPage = event => {
    this.setState({
      perPage: +event.target.value
    });
  };

  getPage = p => {
    const { data, perPage } = this.state;
    if (!data) return [];
    return data.slice(p * perPage, p * perPage + perPage);
  };

  getSelected = () => {
    const {
      name,
      formik: { values }
    } = this.props;
    return getIn(values, name) || [];
  };

  putSelected = k => {
    const {
      name,
      formik: { setFieldTouched, setFieldValue }
    } = this.props;
    const selected = this.getSelected();
    setFieldTouched(name, true);
    let newSelected;
    if (selected.includes(k)) {
      newSelected = selected.filter(v => v !== k);
    } else {
      newSelected = [...selected, k];
    }
    setFieldValue(name, newSelected);
  };

  putAll = () => {
    const { data } = this.state;
    if (!data) return;
    let { keyField } = this.props;
    keyField = keyField || 'id';
    const selected = this.getSelected();
    const selectedAll =
      selected.length > 0 && selected.length === (data || []).length;
    const {
      name,
      formik: { setFieldTouched, setFieldValue }
    } = this.props;
    setFieldTouched(name, true);
    if (selectedAll) {
      setFieldValue(name, []);
    } else {
      setFieldValue(
        name,
        // $FlowFixMe
        data.map(d => d[keyField])
      );
    }
  };

  getError() {
    const {
      name,
      formik: { errors }
    } = this.props;
    return getIn(errors, name);
  }

  hasError() {
    const {
      name,
      formik: { touched }
    } = this.props;
    return getIn(touched, name, false) && this.getError();
  }

  render() {
    const { classes, columns, required, label, helperText } = this.props;

    let { keyField, size } = this.props;

    keyField = keyField || 'id';
    size = size || 'small';

    const {
      data: allData,
      currentPage,
      perPage: rowsPerPage,
      total: totalRows,
      fetching
    } = this.state;

    const isLoading = fetching || !allData;
    const data = this.getPage(currentPage || 0);
    const selected = this.getSelected();
    const selectedAll =
      !isLoading &&
      selected.length > 0 &&
      selected.length === (allData || []).length;
    const selectedAny =
      !isLoading &&
      selected.length > 0 &&
      selected.length < (allData || []).length;
    const finalHelperText = this.hasError() ? this.getError() : helperText;

    return (
      <FormControl
        fullWidth
        required={required}
        className={classes.root}
        error={this.hasError()}
      >
        <FormLabel className={classes.label}>{label}</FormLabel>
        <TableContainer className={classes.container}>
          {isLoading && (
            <div className={classes.loading}>
              <LinearProgress />
            </div>
          )}
          <Table stickyHeader size={size}>
            <TableHeader
              hasCheckbox
              handleSelect={this.putAll}
              selectedAll={selectedAll}
              selectedAny={selectedAny}
              columns={columns}
              sorting={{}}
              sortable={false}
              changeSorting={() => undefined}
            />
            <TableBody
              hasCheckbox
              handleSelect={this.putSelected}
              selectedItems={selected}
              data={data}
              keyField={keyField}
              columns={columns}
              actions={[]}
              size={size}
            />
          </Table>
        </TableContainer>
        <TablePagination
          rowsPerPageOptions={[10, 20, 50, 100]}
          component="div"
          count={totalRows || 0}
          rowsPerPage={rowsPerPage || 10}
          page={currentPage || 0}
          onChangePage={this.handleChangePage}
          onChangeRowsPerPage={this.handleChangeRowsPerPage}
        />
        {!!finalHelperText && (
          <FormHelperText>{finalHelperText}</FormHelperText>
        )}
      </FormControl>
    );
  }
}

export default connect(withStyles(styles)(TableField));
