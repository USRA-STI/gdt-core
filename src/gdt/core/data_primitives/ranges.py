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

__all__ = ['Range', 'TimeRange', 'EnergyRange']

class Range():
    """A primitive class defining a range
    
    Parameters:
        low (float): The low end of the range
        high (float): The high end of the range
    """
    def __init__(self, low, high):
        if high >= low:
            self._low = low
            self._high = high
        else:
            self._low = high
            self._high = low

    @property
    def center(self):
        """(float): The center of the range"""
        return (self._high + self._low) / 2.0

    @property
    def width(self):
        """(float): The width of the range"""
        return self._high - self._low

    def as_tuple(self):
        """Return the range as a tuple.
        
        Returns:  
            (float, float)
        """
        return (self._low, self._high)

    def contains(self, value, inclusive=True):
        """Determine if the range contains a value.
        
        Args:
            value (float): The input value to check
            inclusive (bool, optional): 
                If True, then includes the edges of the range for the check, 
                otherwise it is edge-exclusive. Default is True.
        
        Returns:           
            bool: True if the value is in the range, False otherwise
        """
        if inclusive:
            test = (value <= self._high) and (value >= self._low)
        else:
            test = (value < self._high) and (value > self._low)

        if test:
            return True
        else:
            return False

    @classmethod
    def intersection(cls, range1, range2):
        """Return a new Range that is the intersection of two input Ranges.  
        If the input Ranges do not intersect, then None is returned.

        Args:
            range1 (:class:`Range`): A range used to calculate the intersection
            range2 (:class:`Range`): Another range used to calculate the 
                                     intersection
        
        Returns:           
            :class:`Range`: The intersected range
        """
        # test if one low is inside the other range
        if range1.contains(range2._low):
            lower = range1
            upper = range2
        elif range2.contains(range1._low):
            lower = range2
            upper = range1
        else:
            return None

        # do the merge
        low = upper._low
        high = np.min((lower._high, upper._high))
        obj = cls(low, high)
        return obj

    @classmethod
    def union(cls, range1, range2):
        """Return a new Range that is the union of two input Ranges
        
        Args:
            range1 (:class:`Range`): A range used to calculate the union
            range2 (:class:`Range`): Another range used to calculate the union
        
        Returns:           
            :class:`Range`: The unionized range
        """
        low = np.min((range1._low, range2._low))
        high = np.max((range1._high, range2._high))
        obj = cls(low, high)
        return obj

    def translate(self, value):
        """Returns a new Range that is the translated (shifted) by the given value

        Args:
            value (float): The value to shift the range by

        Returns:
            :class:`Range`: The shifted range
        """
        return self.__class__(self._low + value, self._high + value)

    def __eq__(self, other):
        return (self._low == other._low) and (self._high == other._high)

    def __repr__(self):
        return '<{0}: ({1}, {2})>'.format(self.__class__.__name__, self._low,
                                         self._high)

    def _assert_float(self, value):
        try:
            value = float(value)
        except:
            raise TypeError('value must be a float')
        return value

# mark: TODO: add units to TimeRange/EnergyRange
class TimeRange(Range):
    """A primitive class defining a time range
    
    Parameters:
        tstart (float): The start time of the range
        tstop (float): The end time of the range
    """
    def __init__(self, tstart, tstop):
        tstart = self._assert_float(tstart)
        tstop = self._assert_float(tstop)
        super().__init__(tstart, tstop)

    @property
    def duration(self):
        """(float): The duration of the time range. 
        Alias for ``width``."""
        return self.width

    @property
    def tstart(self):
        """(float): The start time of the range"""
        return self._low

    @property
    def tstop(self):
        """(float): The end time of the range"""
        return self._high


class EnergyRange(Range):
    """A primitive class defining an energy range
    
    Parameters:
        emin (float): The low end of the energy range
        emax (float): The high end of the energy range
    """
    def __init__(self, emin, emax):
        emin = self._assert_float(emin)
        emax = self._assert_float(emax)
        super().__init__(emin, emax)

    @property
    def emax(self):
        """(float): The high end of the energy range"""
        return self._high

    @property
    def emin(self):
        """(float): The low end of the energy range"""
        return self._low

    @property
    def log_center(self):
        """log_center (float): The logarithmic center of the energy range"""
        return np.sqrt(self.emin * self.emax)
