from abc import abstractmethod, ABC
from typing import Any, List, Union, Dict, Callable
import copy

import numpy as np
from astropy.io.fits import FITS_rec

from astropy.time import Time
from astropy.table import Table, Row
from astropy.timeseries import TimeSeries
from astropy.utils.data_info import MixinInfo
from scipy.interpolate import interp1d
from scipy.spatial.transform import Slerp

TableSpec = Union[np.ndarray, Dict[str, Any], Table, FITS_rec]
Numbers = Union[float, int, List[float], List[int], tuple, np.ndarray]


def assert_array_valid(value: np.ndarray, rowsize: int, name: str = 'value', dtype: object = None):
    """Validate the array by checking its size and optionally its dtype.

    Will raise exception if the array doesn't satisfy a condition.
    """
    # Validate array shape
    if rowsize is None or rowsize == 0:
        raise ValueError('The row size specified can not be None or 0.')
    if (value.ndim < 2 and value.size != rowsize != 1) or (value.ndim == 2 and value.shape[1] != rowsize):
        raise ValueError(f'{name} must be an array with a shape of ({rowsize},) or (n, {rowsize}).')
    if value.ndim > 2:
        raise ValueError(f'{name} must not have more than 2 dimensions.')

    # Validate array type
    if dtype is not None and not np.issubdtype(value.dtype, dtype):
        raise ValueError(f'{name} is required to be of type {dtype}.')


def as_times(var: Union[Time, TimeSeries]):
    if isinstance(var, TimeSeries):
        return var['time']
    elif isinstance(var, list) and isinstance(var[0], Time):
        return Time(var)
    elif isinstance(var, Time) and var.isscalar:
        return Time([var])
    elif not isinstance(var, Time):
        raise TypeError('Input value is not a Time object, TimeSeries, or a collection of Time objects.')
    return var


class ArrayBase():
    info = MixinInfo()
    _rowsize: int = None
    _prefix: str = ''
    _seperator: str = ', '
    _dtype: object = np.float64
    _name: str = 'value'

    def __init__(self, value: Any):
        self._orig_shape = None
        self._value = value

    @property
    def _value(self) -> np.ndarray:
        """(np.array): The value in its original shape."""
        return self._array.reshape(self._orig_shape)

    @_value.setter
    def _value(self, value: Any):
        value = np.asarray(value, dtype=self._dtype)
        assert_array_valid(value, self._rowsize, name=self._name, dtype=self._dtype)
        self._orig_shape = value.shape
        self._array = value.reshape(value.size // self._rowsize, self._rowsize)

    def _obj(self):
        """Create a new instance of the class with any info attributes transferred"""
        obj = super().__new__(self.__class__)
        obj._orig_shape = self._orig_shape
        if 'info' in self.__dict__:  # pragma: no cover
            # Coverage doesn't seem to pick this up despite my manually making sure it actually runs in the debugger.
            obj.info = self.info
        return obj

    def copy(self):
        """Returns a copy of the object.
        For compatibility with Astropy frames."""
        obj = copy.copy(self)
        return obj

    @property
    def isscalar(self):
        """(bool): Is or is not scalar"""
        return self.shape == ()

    @property
    def shape(self):
        """(tuple): The shape of the array"""
        if (self._rowsize == 1 and self._orig_shape == ()) or \
                (self._rowsize > 1 and self._orig_shape == (self._rowsize,)):
            return ()
        return self._value.shape[0],

    def __len__(self):
        if self.isscalar:
            raise TypeError(f'object of type {self.__class__.__name__} containing a scalar has no len().')
        else:
            return self._value.shape[0]

    def __getitem__(self, item):
        if self.isscalar:
            raise TypeError(f'object of type {self.__class__.__name__} containing a scalar is not subscriptable.')
        else:
            new = self._obj()
            new._value = self._value[item]
        return new

    def __array__(self, dtype: object = None):
        return np.asarray(self._value, dtype=dtype)

    def _scale(self, value: np.ndarray):
        """Scales the resulting singular value based on the original shape of the internal array.

        Note:
        This only applies to 2-d results that have a single row (1, N) or 1-d results that have a single element (1,).
        Otherwise, the value is returned unchanged
        """
        if self.isscalar:
            if value.ndim == 1 and value.shape == (1,):
                value = value[0]
            elif value.ndim == 2 and value.shape[0] == 1:
                value = value.reshape(value.shape[1])
        return value

    def __str__(self):
        return f'{self._prefix}{np.array2string(self._value, prefix=self._prefix, separator=self._seperator)}'

    def __repr__(self):  # pragma: no coverage
        prefix = f'<{self.__class__.__name__} {self._prefix} '
        return f'{prefix}{np.array2string(self._value, prefix=prefix, separator=self._seperator)} >'


class InterpolationMixin(ABC):
    _interp: Union[interp1d, Slerp, None]
    info: MixinInfo

    def _get_times(self, times) -> Time:
        if times is None:
            try:
                times = self.info.parent_table['time']
            except (AttributeError, KeyError, TypeError):
                raise ValueError('Parent table values are not available. You may be trying to perform operation'
                                 ' with a temporary slice and need to assign the rows in a variable, or as a'
                                 ' stand alone data object without specifying obstimes.')
        return as_times(times)

    def init_interp(self, obstimes: Time):
        """Initialize the interpolation function.

        Args:
            obstimes: Times for the observation data (aligned with the data).
        """
        obstimes = self._get_times(obstimes)
        self._interp = self._setup_interpolation(obstimes.unix_tai)

    @abstractmethod
    def _setup_interpolation(self, unix_tai: Union[float, np.ndarray]) -> Union[interp1d, Slerp]:  # pragma: no cover
        pass

    def at(self, times: Union[Time, TimeSeries], *, obstimes: Union[Time, TimeSeries] = None) -> np.ndarray:
        """Retrieve the interpolated value for the specified times.

        Args:
            times : desired Time(s) in astropy.Time
            obstimes : Times of each ECI location
        Returns:
            np.ndarray aligned with the given times
        """
        if self._interp is None:
            self.init_interp(obstimes)

        return self._interp(as_times(times).unix_tai)


class Interval():
    info = MixinInfo()

    def __init__(self, start: Any, stop: Any):
        self._table = Table(data={
            'start': as_times(start),
            'stop': as_times(stop)
        })
        self.last_inclusive = False

    @classmethod
    def from_table(cls, table: Table,
                   selection: Union[np.ndarray, List[bool], Callable[[Row], bool]] = None) -> 'Interval':
        """Returns an interval containing the earliest and latest time values for the given condition."""
        if selection is None:
            return cls(table['time'].min(), table['time'].max())
        start = []
        stop = []
        if callable(selection):
            previous_state = selection(table[0])
        else:
            previous_state = selection[0]
        if previous_state:
            start.append(table[0]['time'])
        row_index = 1
        last_inclusive = False
        for row in table[1:]:
            if callable(selection):
                state = selection(row)
            else:
                state = selection[row_index]
            if state and not previous_state:
                start.append(row['time'])
            elif not state and previous_state:
                stop.append(row['time'])
            previous_state = state
            row_index += 1
        if len(start) > len(stop):
            stop.append(table[-1]['time'])
            last_inclusive = True
        obj = cls(start, stop)
        obj.last_inclusive = last_inclusive
        return obj

    def _obj(self) -> 'Interval':
        """Create a new instance of the class with any info attributes transferred"""
        obj = super().__new__(self.__class__)
        if 'info' in self.__dict__:  # pragma: no cover
            # Coverage doesn't seem to pick this up despite my manually making sure it actually runs in the debugger.
            obj.info = self.info
        return obj

    def __len__(self):
        return len(self._table)

    @property
    def shape(self) -> tuple:
        return len(self._table),

    @property
    def isscalar(self) -> bool:
        """This is here for compatibility with some of Astropy's functions. It always return False."""
        return False

    def __getitem__(self, item):
        new = self._obj()
        new._table = self._table[item]
        return new

    @property
    def length(self) -> Any:
        return self.stop - self.start

    @property
    def start(self) -> Time:
        return self._table['start']

    @property
    def stop(self) -> Time:
        return self._table['stop']

    @property
    def center(self) -> Any:
        delta = self.length / 2.0
        return self.start + delta

    def in_interval(self, time: Time) -> Union[bool, np.ndarray]:
        for a, b in zip(self.start[:-1], self.stop[:-1]):
            if a <= time < b:
                return True
        if self.last_inclusive:
            if self.start[-1] <= time <= self.stop[-1]:
                return True
        else:
            if self.start[-1] <= time < self.stop[-1]:
                return True

        return False

    def __contains__(self, item):
        return self.in_interval(item)

    def __str__(self):
        lines = []
        for a, b in zip(self.start, self.stop):
            lines.append(f'({a}, {b})')
        return '[' + ',\n '.join(lines) + ']'

    def _repr_html_(self):  # pragma: no cover
        return self._table._repr_html_()


def set_by_dict(obj, values: Dict[str, Any], no_private: bool = False):
    """Sets the object attributes using a dictionary"""
    if no_private:
        valid_args = [p for p in dir(obj) if not p.startswith('_')]
    else:
        valid_args = dir(obj)

    # do the actual setting.
    for key, val in values.items():
        if key in valid_args:
            setattr(obj, key, val)
        else:
            raise AttributeError("{} is not a valid attribute".format(key))
