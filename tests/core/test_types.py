#  CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
#  Contract No.: CA 80MSFC17M0022
#  Contractor Name: Universities Space Research Association
#  Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
#  Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
#  Developed by: William Cleveland and Adam Goldstein
#                Universities Space Research Association
#                Science and Technology Institute
#                https://sti.usra.edu
#
#  Developed by: Daniel Kocevski
#                National Aeronautics and Space Administration (NASA)
#                Marshall Space Flight Center
#                Astrophysics Branch (ST-12)
#
#   Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
#   in compliance with the License. You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software distributed under the License
#  is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
#  implied. See the License for the specific language governing permissions and limitations under the
#  License.
#
from typing import Union, Any

import numpy as np
import pytest
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.timeseries import TimeSeries
from scipy.interpolate import interp1d
from scipy.spatial.transform import Slerp

from gdt.core.types import ArrayBase, assert_array_valid, as_times, InterpolationMixin, Interval, set_by_dict


class TestAssertValidArray:
    @pytest.mark.parametrize("test_array, rowsize, dtype", [
        (1, 1, int),
        ([1], 1, int),
        ([[1], [2], [3]], 1, int),
        ([1., 2., 3.], 3, float),
        ([[1., 2., 3.], [4., 5., 6.]], 3, float)
    ])
    def test_assert_array_valid(self, test_array, rowsize, dtype):
        a = np.asarray(test_array)
        assert_array_valid(a, rowsize=rowsize, dtype=dtype)

    def test_assert_array_invalid_incorrect_rowsize(self):
        a = np.asarray([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(ValueError, match=r'.*shape of \(4,\) or \(n, 4\).*'):
            assert_array_valid(a, rowsize=4, dtype=int)

    def test_assert_array_invalid_more_than_2d(self):
        a = np.asarray([[[1, 2, 3], [4, 5, 6]]])
        with pytest.raises(ValueError, match=r'.*must not have more than 2 dimensions.*'):
            assert_array_valid(a, rowsize=4, dtype=int)

    def test_assert_array_invalid_dtype(self):
        a = np.asarray([['A', 'B', 'C'], [4, 5, 6]])
        with pytest.raises(ValueError, match=r'.*value is required to be of type.*'):
            assert_array_valid(a, rowsize=3, dtype=int)


class TestAsTimes:
    @pytest.mark.parametrize("in_value, expected_len, expected_value", [
        (Time('2008-06-11T16:05:00', format='isot', scale='utc'), 1,
         Time(['2008-06-11T16:05:00'], format='isot', scale='utc')),
        ([Time('2008-06-11T16:05:00', format='isot', scale='utc'),
          Time('2008-06-11T16:05:01', format='isot', scale='utc'),
          Time('2008-06-11T16:05:02', format='isot', scale='utc')], 3,
         Time(['2008-06-11T16:05:00', '2008-06-11T16:05:01', '2008-06-11T16:05:02'], format='isot', scale='utc')),
        (
                TimeSeries(time_start=Time('2008-06-11T16:05:00', format='isot', scale='utc'), time_delta=1 * u.s,
                           n_samples=3),
                3,
                Time(['2008-06-11T16:05:00', '2008-06-11T16:05:01', '2008-06-11T16:05:02'], format='isot',
                     scale='utc'))])
    def test_as_times(self, in_value, expected_len, expected_value):
        result = as_times(in_value)
        assert isinstance(result, Time)
        assert not result.isscalar
        assert len(result) == expected_len
        assert np.all(result.utc.isot == expected_value.utc.isot)

    def test_invalid_time(self):
        with pytest.raises(TypeError, match=r'.*is not a Time object.*'):
            as_times('Bad Value')


class TestArrayBase:
    def test_instantiation_base(self):
        with pytest.raises(ValueError, match=r"The row size specified can not be None or 0."):
            ArrayBase(1)

    def test_instantiation_with_none(self):
        class Data(ArrayBase):
            _rowsize = 4

        with pytest.raises(ValueError, match=r'.*\(4,\) or \(n, 4\).*'):
            Data(None)

    def test_instantiation_with_wrong_shape(self):
        class Data(ArrayBase):
            _rowsize = 4

        with pytest.raises(ValueError, match=r'.*\(4,\) or \(n, 4\).*'):
            Data([1, 2, 3])

        with pytest.raises(ValueError, match=r'.*\(4,\) or \(n, 4\).*'):
            Data([1, 2, 3, 4, 5])

        with pytest.raises(ValueError, match=r'.*inhomogeneous shape.*'):
            Data([[1, 2, 3, 4], [1, 2, 3]])

    def test_instantiation_row_of_one_single_value_0d(self):
        class Data(ArrayBase):
            _rowsize = 1

        d = Data(128.125)
        assert d.isscalar
        assert d._array.shape == (1, 1)
        assert np.asarray(d).shape == ()
        assert d.shape == ()
        with pytest.raises(TypeError):
            len(d)
        with pytest.raises(TypeError):
            d[0]

    def test_instantiation_row_of_one_single_value_1d(self):
        class Data(ArrayBase):
            _rowsize = 1

        d = Data([128.125])
        assert not d.isscalar
        assert d._array.shape == (1, 1)
        assert np.asarray(d).shape == (1,)
        assert d.shape == (1,)
        assert len(d) == 1

    def test_instantiation_multiple_values_of_one(self):
        class Data(ArrayBase):
            _rowsize = 1

        d = Data([1, 2, 3, 4, 5])
        assert not d.isscalar
        assert d._array.shape == (5, 1)
        assert np.asarray(d).shape == (5,)
        assert d.shape == (5,)
        assert len(d) == 5

    def test_instantiation_multiple_values_of_one_2d(self):
        class Data(ArrayBase):
            _rowsize = 1

        d = Data([[1], [2], [3], [4], [5]])
        assert not d.isscalar
        assert d._array.shape == (5, 1)
        assert np.asarray(d).shape == (5, 1)
        assert d.shape == (5,)
        assert len(d) == 5

    def test_instantiation_of_four_as_1d(self):
        class Data(ArrayBase):
            _rowsize = 4

        with pytest.raises(ValueError):
            Data([1, 2, 3])

        with pytest.raises(ValueError):
            Data([1, 2, 3, 4, 5])

        d = Data([1, 2, 3, 4])
        assert d.isscalar
        assert d._array.shape == (1, 4)
        assert np.asarray(d).shape == (4,)
        assert d.shape == ()
        with pytest.raises(TypeError):
            len(d)
        with pytest.raises(TypeError):
            d[0]

    def test_instantiation_of_multiple_values_of_four(self):
        class Data(ArrayBase):
            _rowsize = 4

        with pytest.raises(ValueError):
            Data([[1, 2, 3, 4], [1, 2, 3]])

        with pytest.raises(ValueError):
            Data([[1, 2, 3, 4, 5], [1, 2, 3, 4]])

        d = Data([[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]])
        assert not d.isscalar
        assert d._array.shape == (3, 4)
        assert np.asarray(d).shape == (3, 4)
        assert d.shape == (3,)
        assert len(d) == 3

    def test_retrieving_items_of_one(self):
        class Data(ArrayBase):
            _rowsize = 1

        d = Data([1, 2, 3, 4, 5])
        assert len(d) == 5
        assert d.shape == (5,)

        e = d[1]
        assert isinstance(e, Data)
        assert e.isscalar
        assert np.asarray(e) == np.array(2)

        e = d[:3]
        assert isinstance(e, Data)
        assert not e.isscalar
        assert np.all(np.asarray(e) == np.array([1, 2, 3]))
        assert np.asarray(e).shape == (3,)

    def test_retrieving_items_of_one_2d(self):
        class Data(ArrayBase):
            _prefix = '(x)'
            _rowsize = 1

        d = Data([[1], [2], [3], [4], [5]])
        assert len(d) == 5
        assert d.shape == (5,)
        assert np.asarray(d).shape == (5, 1)
        assert str(d) == '(x)[[1.],\n    [2.],\n    [3.],\n    [4.],\n    [5.]]'

        e = d[1]
        assert isinstance(e, Data)
        assert not e.isscalar
        assert np.asarray(e) == np.array([2])
        assert str(e) == '(x)[2.]'

        e = d[:3]
        assert isinstance(e, Data)
        assert not e.isscalar
        assert np.all(np.asarray(e) == np.array([[1], [2], [3]]))
        assert np.asarray(e).shape == (3, 1)

    def test_retrieving_items_of_four(self):
        class Data(ArrayBase):
            _rowsize = 4

        d = Data([[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6], [4, 5, 6, 7]])
        assert len(d) == 4
        assert d.shape == (4,)

        e = d[1]
        assert isinstance(e, Data)
        assert e.isscalar
        assert np.all(np.asarray(e) == np.array([2, 3, 4, 5]))
        assert np.asarray(e).shape == (4,)
        assert e.shape == ()
        with pytest.raises(TypeError):
            len(e)

        e = d[1:]
        assert isinstance(e, Data)
        assert not e.isscalar
        assert np.all(np.asarray(e) == np.array([[2, 3, 4, 5], [3, 4, 5, 6], [4, 5, 6, 7]]))
        assert np.asarray(e).shape == (3, 4)
        assert e.shape == (3,)
        assert len(e) == 3

    def test_scale_method(self):
        class Data(ArrayBase):
            _rowsize = 4

            @property
            def val_1d(self):
                return self._scale(self._array[:, 2])

            @property
            def val_2d(self):
                return self._scale(self._array[:, :2])

        d = Data([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        assert not d.isscalar
        assert d._array.shape == (3, 4)
        assert np.asarray(d).shape == (3, 4)
        assert d.shape == (3,)
        assert len(d) == 3
        assert np.all(d.val_1d == [3, 7, 11])
        assert np.all(d.val_2d == [[1, 2], [5, 6], [9, 10]])

        e = d[1]
        assert e.isscalar
        assert e._array.shape == (1, 4)
        assert np.asarray(e).shape == (4,)
        assert e.shape == ()
        assert e.val_1d == 7
        assert np.all(e.val_2d == [5, 6])


class TestInterpolationMixin:
    class DataClass(InterpolationMixin, ArrayBase):
        _rowsize = 1

        def __init__(self, value: Any):
            super().__init__(value)
            self._interp = None

        def _setup_interpolation(self, unix_tai: Union[float, np.ndarray]) -> Union[interp1d, Slerp]:
            return interp1d(unix_tai, self._value)

    def test_interpolation(self):
        obstime = Time(['2008-06-11T16:05:00', '2008-06-11T16:05:01', '2008-06-11T16:05:02', '2008-06-11T16:05:03',
                        '2008-06-11T16:05:09', '2008-06-11T16:05:10'], format='isot', scale='utc')

        d = self.DataClass((11, 11.5, 12, 12.5, 15.5, 16))
        ts = TimeSeries(time_start=Time('2008-06-11T16:05:04', format='isot', scale='utc'), time_delta=1 * u.s,
                        n_samples=5)

        vals = d.at(ts, obstimes=obstime)
        assert np.all(vals == (13, 13.5, 14, 14.5, 15))

    def test_interpolation_within_table(self):
        table = TimeSeries(data={
            'data': self.DataClass((11, 11.5, 12, 12.5, 15.5, 16))
        }, time=Time(['2008-06-11T16:05:00', '2008-06-11T16:05:01', '2008-06-11T16:05:02', '2008-06-11T16:05:03',
                      '2008-06-11T16:05:09', '2008-06-11T16:05:10'], format='isot', scale='utc'))

        ts = TimeSeries(time_start=Time('2008-06-11T16:05:04', format='isot', scale='utc'), time_delta=1 * u.s,
                        n_samples=5)

        vals = table['data'].at(ts)
        assert np.all(vals == (13, 13.5, 14, 14.5, 15))

    def test_interpolation_without_table_or_time(self):
        d = self.DataClass((11, 11.5, 12, 12.5, 15.5, 16))
        ts = TimeSeries(time_start=Time('2008-06-11T16:05:04', format='isot', scale='utc'), time_delta=1 * u.s,
                        n_samples=5)

        with pytest.raises(ValueError):
            d.at(ts)


@pytest.mark.parametrize("data, num_intervals, in_time, out_time", [
    ((True, True, True, True, False, False, False, False, True, True), 2, (0, 1, 2, 3, 8, 9), (4, 5, 6, 7)),
    ((False, False, True, True, False, False, True, True, False, False), 2, (2, 3, 6, 7), (0, 1, 4, 5, 8, 9)),
    ((False, True, False, True, False, True, False, True, False, True), 5, (1, 3, 5, 7, 9), (0, 2, 4, 6, 8)),
    ((True, False, True, False, True, False, True, False, True, False), 5, (0, 2, 4, 6, 8), (1, 3, 5, 7, 9))
])
class TestInterval:

    def test_interval_from_table_lambda(self, data, num_intervals, in_time, out_time):
        t0 = Time('2022-01-01T00:00:00', format='isot', scale='utc')
        ts = TimeSeries(data={'ok': data}, time_start=t0, time_delta=1 * u.s, n_samples=len(data))

        interval = Interval.from_table(ts, lambda x: x['ok'] == True)
        assert len(interval) == num_intervals
        assert interval.shape == (num_intervals,)
        # Check all the times inside intervals
        for x in in_time:
            t = t0 + TimeDelta(x, format='sec')
            assert t in interval
        # Check all the times outside intervals
        for x in out_time:
            t = t0 + TimeDelta(x, format='sec')
            assert t not in interval

    def test_interval_from_table_mask(self, data, num_intervals, in_time, out_time):
        t0 = Time('2022-01-01T00:00:00', format='isot', scale='utc')
        ts = TimeSeries(data={'ok': data}, time_start=t0, time_delta=1 * u.s, n_samples=len(data))

        interval = Interval.from_table(ts, ts['ok'] == True)
        assert len(interval) == num_intervals
        # Check all the times inside intervals
        for x in in_time:
            t = t0 + TimeDelta(x, format='sec')
            assert t in interval
        # Check all the times outside intervals
        for x in out_time:
            t = t0 + TimeDelta(x, format='sec')
            assert t not in interval

    def test_interval_no_selection(self, data, num_intervals, in_time, out_time):
        t0 = Time('2022-01-01T00:00:00', format='isot', scale='utc')
        te = t0 + TimeDelta(len(data) - 1, format='sec')
        ts = TimeSeries(data={'ok': data}, time_start=t0, time_delta=1 * u.s, n_samples=len(data))

        interval = Interval.from_table(ts)
        assert len(interval) == 1

        assert not interval.isscalar
        assert interval.start == t0
        assert interval.stop == te
        assert interval.center == t0 + (te - t0) / 2
        assert np.allclose(interval.length.sec, [9.])
        assert interval[0].start == t0
        assert interval[0].stop == te

        assert str(interval) == f'[({t0.utc.isot}, {te.utc.isot})]'


class TestSetByDict:
    class MockClass:
        a = 12

        def __init__(self):
            self.b = 0
            self._c = None

        @property
        def c(self):
            return self._c

        @c.setter
        def c(self, val):
            self._c = val

    def test_set_values(self):
        mock = self.MockClass()
        set_by_dict(mock, {'a': 14, 'b': 20, 'c': 'Hello'})
        assert mock.a == 14
        assert mock.b == 20
        assert mock.c == 'Hello'

        set_by_dict(mock, {'_c': 40})
        assert mock.c == 40

    def test_set_values_no_private(self):
        mock = self.MockClass()
        set_by_dict(mock, {'a': 14, 'b': 20, 'c': 'Hello'}, no_private=True)
        assert mock.a == 14
        assert mock.b == 20
        assert mock.c == 'Hello'

        with pytest.raises(AttributeError):
            set_by_dict(mock, {'_c': 40}, no_private=True)

    def test_set_value_not_an_attribute(self):
        mock = self.MockClass()
        with pytest.raises(AttributeError):
            set_by_dict(mock, {'fake': 100})
