# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract Nos.: CA 80MSFC17M0022 / 80NSSC24M0035
# Contractor Name: Universities Space Research Association
# Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
# Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland and Adam Goldstein
#               Universities Space Research Association
#               Science and Technology Institute
#               https://sti.usra.edu
#
# Developed by: Daniel Kocevski
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#
import numpy as np
from .ranges import Range, TimeRange, EnergyRange

__all__ = ['Intervals', 'Gti', 'Ebounds']


class Intervals():
    """A primitive class defining a set of intervals or ranges.
    
    An interval may be accessed by using indices.
    
    Parameters:
        interval (:class:`Range`, optional): An interval to initialize with
    """
    _range_class = Range
    def __init__(self, interval=None):
        self._intervals = None
        if interval is not None:
            self._intervals = [interval]

    @property
    def intervals(self):
        """(list): The list of intervals"""
        return self._intervals

    @property
    def num_intervals(self):
        """(int): The number of intervals"""
        return len(self._intervals)

    @property
    def range(self):
        """(float, float): The full range spanned by the intervals"""
        if self.num_intervals > 0:
            return (self._intervals[0].as_tuple()[0],
                    self._intervals[-1].as_tuple()[1])

    def as_list(self):
        """Return the intervals as a list of tuples.
        
        Returns:     
            [(float, float), ...]
        """
        interval_list = [interval.as_tuple() for interval in self._intervals]
        return interval_list

    def contains(self, value, inclusive=True):
        """Determine if the intervals contains a value.
        
        Args:
            value (float): The input value to check
            inclusive (bool, optional): 
                If True, then includes the edges of the range for the check, 
                otherwise it is edge-exclusive. Default is True.
        
        Returns:           
            bool: True if the value is in the intervals, False otherwise
        """
        test = [interval.contains(value, inclusive=inclusive) for interval \
                in self._intervals]
        return any(test)

    def high_edges(self):
        """Return a list of the high edges.
        
        Returns:     
            (list)
        """
        edges = [interval._high for interval in self._intervals]
        return edges

    def index(self, value):
        """Return the index of the interval that contains the value.
        
        Note :
            If the value is precisely on the boundary of two intervals, the
            first interval index will be returned.
        
        Args:
            value (float): The input value
        
        Returns:
            (int)
        """
        for i, interval in enumerate(self._intervals):
            if interval.contains(value, inclusive=True):
                return i

    def insert(self, interval):
        """Insert a new interval
        
        Args:
            interval (:class:`Range`): The interval to insert
        """

        # if interval is a duplicate, then skip
        for _interval in self._intervals:
            if interval == _interval:
                return

        # where the new range should be inserted
        idx = [i for i, j in enumerate(self._intervals) if
               j._low <= interval._low]

        # determine if there is overlap with the lower bounding range, and if 
        # so, then merge
        if len(idx) != 0:
            idx = idx[-1]
            if self._intervals[idx].contains(interval._low, inclusive=False):
                the_range = self._intervals.pop(idx)
                interval = type(interval).union(the_range, interval)
            else:
                idx += 1
        else:
            idx = 0

        # determine if there is overlap with the upper bounding range, and if 
        # so, then merge
        if idx < len(self._intervals):
            if self._intervals[idx].contains(interval._high):
                the_range = self._intervals.pop(idx)
                interval = type(interval).union(the_range, interval)

        self._intervals.insert(idx, interval)

    def low_edges(self):
        """Return a list of the low edges.
        
        Returns:     
            (list)
        """
        edges = [interval._low for interval in self._intervals]
        return edges

    @classmethod
    def from_bounds(cls, low_bounds, high_bounds):
        """Create a new Intervals object from a list of lower and upper bounds.
        
        Args:
            low_bounds (list): The lower bounds of the intervals
            high_bounds (list): The upper bounds of the intervals
        
        Returns:           
            :class:`Intervals`
        """
        low_bounds = np.asarray(low_bounds).flatten()
        high_bounds = np.asarray(high_bounds).flatten()
        num = low_bounds.size
        if num != high_bounds.size:
            raise ValueError('low_bounds and high_bounds must be of same size')

        intervals = [cls._range_class(low_bounds[i], high_bounds[i]) \
                     for i in range(num)]
        obj = cls(intervals[0])
        obj._intervals = intervals
        return obj

    @classmethod
    def from_list(cls, interval_list):
        """Create a new Intervals object from a list of tuples.
        
        Args:
            interval_list ([(float, float), ...]):  A list of interval tuples
        
        Returns:           
            :class:`Intervals`
        """
        intervals = [cls._range_class(*interval) for interval in interval_list]
        obj = cls(intervals[0])
        obj._intervals = intervals
        return obj

    @classmethod
    def intersection(cls, intervals1, intervals2):
        """Return a new Intervals object that is the intersection of two 
        existing Intervals objects.
        
        Args:
            intervals1 (:class:`Intervals`): Intervals to be intersected
            intervals2 (:class:`Intervals`): Intervals to be intersected
           
        Returns:
            :class:`Intervals`
        """
        if intervals1.range[0] <= intervals2.range[0]:
            _intervals1 = intervals1
            _intervals2 = intervals2
        else:
            _intervals1 = intervals2
            _intervals2 = intervals1

        new_intervals = []
        for interval1 in _intervals1.intervals:
            for interval2 in _intervals2.intervals:
                if not interval1.contains(interval2._low) and \
                   not interval1.contains(interval2._high):
                    continue
                new_intervals.append(cls._range_class.intersection(interval1,
                                                                   interval2))
        obj = cls()
        obj._intervals = new_intervals
        return obj

    @classmethod
    def merge(cls, intervals1, intervals2):
        """Return a new Intervals object that is a merge of two existing 
        Intervals objects.
        
        Args:
            intervals1 (:class:`Intervals`): Intervals to be merged
            intervals2 (:class:`Intervals`): Intervals to be merged
           
        Returns:
            (:class:`Intervals`)
        """
        interval_list = intervals1.as_list()
        interval_list.extend(intervals2.as_list())
        interval_list = sorted(interval_list)
        obj = cls(cls._range_class(*interval_list.pop(0)))
        for interval in interval_list:
            obj.insert(cls._range_class(*interval))
        return obj

    def __repr__(self):
        s = '<{0}: {1} intervals; range {2}>'.format(self.__class__.__name__,
                                                     self.num_intervals,
                                                     self.range)
        return s

    def __getitem__(self, index):
        if index > self.num_intervals-1:
            raise KeyError('Requested {0}th interval and only {1} intervals '\
                           'exist'.format(index, self.num_intervals))
        return self._intervals[index]


class Gti(Intervals):
    """A primitive class defining a set of Good Time Intervals (GTIs).
    
    An interval may be accessed by using indices.
    
    Parameters:
        interval (:class:`TimeRange`, optional): An interval to initialize with
    """
    _range_class = TimeRange

    @classmethod
    def from_boolean_mask(cls, times, mask):
        """Create a new GTI object from a list of times and a Boolean mask
        Splits the boolean mask into segments of contiguous values and applies
        to array of times to create a GTI object.
        
        Args:
            times (np.array): An array of times
            mask (np.array(dtype=bool)): The boolean array. Must be the same 
                                         size as times.
        
        Returns:           
            :class:`Gti`: The new GTI object
        """
        times = np.asarray(times)
        mask = np.asarray(mask)

        # split a boolean mask array into segments based on True/False
        indices = np.nonzero(mask[1:] != mask[:-1])[0] + 1
        time_segs = np.split(times, indices)
        mask_segs = np.split(mask, indices)

        # retrieve the start and stop times for the "on" intervals
        segs = []
        numsegs = len(indices) + 1
        for i in range(numsegs):
            if mask_segs[i][0]:
                segs.append((time_segs[i][0], time_segs[i][-1]))

        # if mask is all True or all False
        if len(segs) == 0:
            if mask[0]:
                segs = [(times.min(), times.max())]
            else:
                return None

        return cls.from_list(segs)


class Ebounds(Intervals):
    """A primitive class defining a set of energy bounds.
    
    An interval may be accessed by using indices.
    
    Parameters:
        interval (:class:`EnergyRange`, optional): An interval to initialize with
    """
    _range_class = EnergyRange
