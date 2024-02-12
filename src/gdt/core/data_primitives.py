# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
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
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import copy

__all__ = ['Range', 'TimeRange', 'EnergyRange', 'Intervals', 'Gti', 'Ebounds',
           'EventList', 'Bins', 'ExposureBins', 'ChannelBins', 'TimeBins', 
           'EnergyBins', 'TimeChannelBins', 'TimeEnergyBins', 'ResponseMatrix', 
           'Parameter']

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
    

class EventList():
    """A primitive class defining an event list.
    
    Parameters:
        times (np.array): The array of event times
        channels (np.array): The corresponding array of associated energy channel
        ebounds (:class:`Ebounds`, optional): The energy bounds mapped to 
                                              energy channels
    """
    def __init__(self, times=None, channels=None, ebounds=None):
        
        if (times is not None) and (channels is not None):
            times = np.asarray(times).flatten()
            channels = np.asarray(channels).flatten()
            if times.size != channels.size:
                raise ValueError('times and channels arrays must be same length')
        else:
            times = np.array([], dtype=float)
            channels = np.array([], dtype=int)
        
        events = zip(*(times, channels))
        self._events = np.array(list(events), 
                                dtype=[('TIME', times.dtype.type),
                                       ('PHA', channels.dtype.type)])
        
        if ebounds is not None:
            if not isinstance(ebounds, Ebounds):
                raise TypeError('ebounds must be of type Ebounds')
        self._ebounds = ebounds
        
    @property
    def channel_range(self):
        """(int, int): The range of the channels in the list"""
        if self.size > 0:
            return int(np.min(self.channels)), int(np.max(self.channels))

    @property
    def channels(self):
        """(np.array): The PHA channel array"""
        return self._events['PHA']
    
    @property
    def ebounds(self):
        """(:class:`Ebounds`): The energy bounds of the energy channels.
                               This property can be set.
        """
        return self._ebounds
    
    @ebounds.setter
    def ebounds(self, val):
        if not isinstance(val, Ebounds):
            raise TypeError('ebounds must be an Ebounds object')
        self._ebounds = val
        
    @property
    def emax(self):
        """(float): The maximum energy"""
        if self._ebounds is not None:
            return self._ebounds.range[1]

    @property
    def emin(self):
        """(float): The minimum energy"""
        if self._ebounds is not None:
            return self._ebounds.range[0]

    @property
    def energy_range(self):
        """(float, float): The energy range of the channels"""
        if self.size > 0:
            if self._ebounds is not None:
                chan_min, chan_max = self.channel_range
                return (self._ebounds[chan_min].emin,
                        self._ebounds[chan_max].emax)

    @property
    def num_chans(self):
        """(int): The number of energy channels. 
        Note that not all channels will necessarily have events, especially if a 
        slice is made over energy."""
        if self.ebounds is not None:
            return self._ebounds.num_intervals
        else:
            return (self._events['PHA'].max() - self._events['PHA'].min()) + 1

    @property
    def size(self):
        """(int): The number of events in the list"""
        return self._events.size
    
    @property
    def time_range(self):
        """(float, float): The range of the times in the list"""
        if self.size > 0:
            return self.times.min(), self.times.max()

    @property
    def times(self):
        """(np.array): The event times"""
        return self._events['TIME']
    
    def bin(self, method, *args, tstart=None, tstop=None,
            event_deadtime=2.6e-6, overflow_deadtime=1.0e-5, **kwargs):
        """Bin the EventList in time given a binning function and return a
        2D time-energy channel histogram. If the ebounds energy calibration is
        set, returns a :class:`TimeEnergyBins` object, otherwise returns a
        :class:`TimeChannelBins` object.
        
        The binning function should take as input an array of times as well
        as a tstart and tstop keywords for partial list binning.  Additional 
        arguments and keyword arguments specific to the function are allowed.
        The function should return an array of time edges for the bins, such
        that, for `n` bins, there are `n` + 1 edges.
        
        Args:
            method (<function>): A binning function
            *args: Arguments to be passed to the binning function
            tstart (float, optional): 
                If set, defines the start time of the EventList to be binned, 
                otherwise binning will begin at the time of the first event.
            tstop (float, optional): 
                If set, defines the end time of the EventList to be binned, 
                otherwise binning will end at the time of the last event.
            event_deadtime (float, optional): The deadtime per event in seconds. 
                                              Default is 2.6e-6.
            overflow_deadtime (float, optional): 
                The deadtime per event in the overflow channel in seconds. 
                Default is 1e-5.
            **kwargs: Options to be passed to the binning function
        
        Returns:
            (:class:`TimeEnergyBins` or :class:`TimeChannelBins`)
        """
        if tstart is None:
            tstart = self.time_range[0]
        if tstop is None:
            tstop = self.time_range[1]

        # set the start and stop of the rebinning segment        
        mask = (self.times >= tstart) & (self.times <= tstop)
        bin_times = self.times[mask]
        kwargs['tstart'] = tstart
        kwargs['tstop'] = tstop

        # get the time edges from the binning function and then do the 2d histo
        time_edges = method(bin_times, *args, **kwargs)
        if self.ebounds is not None:
            chan_list = np.arange(self.ebounds.num_intervals + 1)
        else:
            chan_list = np.arange(self.channel_range[1] + 2)
        counts = np.histogram2d(bin_times, self.channels[mask],
                                [time_edges, chan_list])[0]

        # calculate exposure
        lo_edges = time_edges[:-1]
        hi_edges = time_edges[1:]
        overflow_counts = counts[:, -1]
        deadtime = counts.sum(axis=1) * event_deadtime + \
                   overflow_counts * (overflow_deadtime - event_deadtime)
        exposure = (hi_edges - lo_edges) - deadtime
        
        if self.ebounds is None:
            bins = TimeChannelBins(counts, lo_edges, hi_edges, exposure,
                                   chan_list[:-1])
        else:
            bins = TimeEnergyBins(counts, lo_edges, hi_edges, exposure,
                                  self.ebounds.low_edges(), 
                                  self.ebounds.high_edges())
        return bins

    def channel_slice(self, chanlo, chanhi):
        """Perform a slice in energy channels of the EventList and return a 
        new EventList
        
        Args:
            chanlo (int): The start of the channel slice
            chanhi (int): The end of the channel slice
        
        Returns:        
            (:class:`EventList`)
        """
        mask = (self.channels >= chanlo) & (self.channels <= chanhi)
        events = self._events[mask]
        return EventList(times=events['TIME'], channels=events['PHA'],
                         ebounds=copy.deepcopy(self._ebounds))

    def count_spectrum(self, **kwargs):
        """Extract a count spectrum for the EventList or for a segment of the
        EventList. If the ebounds energy calibration is set, returns an 
        :class:`EnergyBins` object, otherwise returns a :class:`ChannelBins`
        object.
        
        Args:
            time_ranges ([(float, float), ...], optional): 
                The time range or time ranges over which to calculate the 
                exposure. If omitted, calculates the total exposure of the data.
            event_deadtime (float, optional): The deadtime per event in seconds. 
                                              Default is 0.
            overflow_deadtime (float, optional): 
                The deadtime per event in the overflow channel in seconds. 
                Default is 0.
        
        Returns:        
            (:class:`EnergyBins` or :class:`ChannelBins`)
        """
        if self.ebounds is not None:
            chan_list = np.arange(self.ebounds.num_intervals + 1)
        else:    
            chan_list = np.arange(self.channel_range[1] + 2)
        counts = np.histogram(self.channels, bins=chan_list)[0]
        exposure = np.full(chan_list.size-1, self.get_exposure(time_ranges=None,
                                                               **kwargs))
        
        if self.ebounds is None:
            bins = ChannelBins.create(counts, chan_list[:-1], exposure)
        else:
            bins = EnergyBins(counts, self.ebounds.low_edges(), 
                              self.ebounds.high_edges(), exposure)
        return bins

    def energy_slice(self, emin, emax):
        """Perform a slice in energy of the EventList and return a new EventList.
        Since energies are binned, an emin or emax falling inside of an energy
        channel bin will include that bin in the slice.
        
        Args:
            emin (float): The start of the energy slice
            emax (float): The end of the energy slice
        
        Returns:        
            (:class:`EventList`)
        """
        if self._ebounds is None:
            raise AttributeError('Ebounds has not been set, therefore the \
                                  the EventList cannot be sliced by energy')
        emins = np.array([ebound.emin for ebound in self._ebounds.intervals])
        emaxs = np.array([ebound.emax for ebound in self._ebounds.intervals])
        mask = (emins < emax) & (emaxs > emin)
        chan_nums = np.arange(self.channel_range[0], self.channel_range[1]+1)
        chan_nums = chan_nums[mask]        
        return self.channel_slice(chan_nums[0], chan_nums[-1])
    
    def get_exposure(self, time_ranges=None, event_deadtime=0.0,
                     overflow_deadtime=0.0):
        """Calculate the total exposure, in seconds, of a time range or time 
        ranges of data.
        
        Args:
            time_ranges ([(float, float), ...], optional): 
                The time range or time ranges over which to calculate the 
                exposure. If omitted, calculates the total exposure of the data.
            event_deadtime (float, optional): The deadtime per event in seconds. 
                                              Default is 0.
            overflow_deadtime (float, optional): 
                The deadtime per event in the overflow channel in seconds. 
                Default is 0.
        
        Returns:
            (float)
        """
        if time_ranges is None:
            time_ranges = [self.time_range]
        try:
            iter(time_ranges[0])
        except:
            time_ranges = [time_ranges]

        exposure = 0.0
        for i in range(len(time_ranges)):
            tstart, tstop = self._assert_range(time_ranges[i])
            tcent = (tstop + tstart) / 2.0
            dt = (tstop - tstart) / 2.0
            mask = (np.abs(self.times - tcent) <= dt)
            # mask = (self.time >= tstart) & (self.time < tstop)
            deadtime = mask.sum() * event_deadtime
            
            if self.ebounds is not None:
                of_chan = self.ebounds.num_intervals-1
            else:
                of_chan = self.channel_range[1]
            
            deadtime += np.sum(self.channels[mask] == of_chan) * \
                        (overflow_deadtime - event_deadtime)
            exposure += (tstop - tstart) - deadtime
        return exposure

    def rebin_energy(self, method, *args, **kwargs):
        """Rebin the energy channels using the specified binning algorithm.  
        This does not change the number of events in the EventList, but changes 
        their assignment to a channel and bins the energy bounds mapping to 
        those channels (if applicable).  A new EventList object is returned.
        
        Args:
            method (<function>): A binning function
            *args: Arguments to be passed to the binning function
            **kwargs: Deadtime arguments passed to get_exposure
        
        Returns:
           (:class:`EventList`)
        """
        # exposure and counts; not really used other than for some specific
        # binning algorithms
        exposure = np.full(self.num_chans, self.get_exposure(**kwargs))
        chans, counts = np.unique(self.channels, return_counts=True)
        if chans.size != self.num_chans:
            counts_fill = np.zeros(self.num_chans)
            counts_fill[chans] = counts
            counts = counts_fill
        edges = np.arange(self.num_chans + 1)
        
        # call the binning algorithm and get the new edges
        _, _, new_edges = method(counts, exposure, edges, *args)

        # re-assign the pha channels based on the new edges
        # and also rebin the ebounds
        new_channels = np.zeros_like(self.channels)
        new_ebounds = []
        for i in range(len(new_edges) - 1):
            emin = new_edges[i]
            emax = new_edges[i + 1]
            mask = (self.channels >= emin) & (self.channels < emax)
            new_channels[mask] = i
            if self.ebounds is not None:
                new_ebounds.append((self.ebounds[emin].emin,
                                    self.ebounds[emax - 1].emax))

        # create the new EventList object with the rebinned channels
        if self.ebounds is not None:
            new_ebounds = Ebounds.from_list(new_ebounds)
        else:
            new_ebounds = None
        obj = EventList(times=np.copy(self.times), channels=new_channels,
                        ebounds=new_ebounds)
        return obj

    def sort_time(self):
        """In-place sort by time.
        """
        idx = np.argsort(self._events['TIME'])
        self._events = self._events[:][idx]

    def sort_channels(self):
        """In-place sort by channel number.
        """
        idx = np.argsort(self._events['PHA'])
        self._events = self._events[:][idx]
    
    def time_slice(self, tstart, tstop):
        """Perform a slice in time of the EventList and return a new EventList
        
        Args:
            tstart (float): The start of the time slice
            tstop (float): The end of the time slice
        
        Returns:        
            (:class:`EventList`)
        """
        mask = (self.times >= tstart) & (self.times <= tstop)
        events = self._events[mask]
        return EventList(times=events['TIME'], channels=events['PHA'],
                         ebounds=copy.deepcopy(self._ebounds))

    @classmethod
    def merge(cls, eventlists, sort=False, force_unique=True):
        """Merge multiple EventLists together in time and optionally sort.
        
        Args:
            eventlist (list of :class:`EventList`): 
                A list containing the EventLists to be merged
            sort (bool, optional): 
                If True, sorts by time after the merge. Default is False.
            force_unique (bool, optional): 
                If True, force all events to be unique via brute force sorting. 
                If False, the EventLists will only be checked and masked for 
                overlapping time ranges. Events can potentially be lost if the 
                merged EventLists contain overlapping times (but not necessarily 
                duplicate events), however this method is much faster.  
                Default is True.
        
        Returns:
            (:class:`EventList`)
        """
                
        # put in time order
        idx = np.argsort([eventlist.times.min() for eventlist in eventlists])
        eventlists = [eventlists[i] for i in idx]
        
        
        # mark: TODO add check for ebounds consistency
        ebounds_idx = 0
        for i in range(len(eventlists)):
            if eventlists[i].ebounds is None:
                continue
            else:
                ebounds_idx = i
                break
        
        new_times = np.array([])
        new_chans = np.array([])
        for eventlist in eventlists:
            # skip empty EventLists
            if eventlist.size == 0:
                continue
            
            # if not forcing to be unique, just make sure there is no time overlap
            if (not force_unique) and (new_times.size > 0):
                mask = (eventlist.times > new_times.max())
                temp_times = eventlist.times[mask]
                temp_chans = eventlist.channels[mask]
            else:
                temp_times = eventlist.times
                temp_chans = eventlist.channels
            new_times = np.concatenate((new_times, temp_times))
            new_chans = np.concatenate((new_chans, temp_chans))

        # force unique: make sure that we only keep unique events (slower)
        if force_unique:
            new_times, uidx = np.unique(new_times, return_index=True)
            new_chans = new_chans[uidx]
            
        obj = cls(times=new_times, channels=new_chans, 
                  ebounds=copy.deepcopy(eventlists[ebounds_idx].ebounds))

        # do a sort
        if sort and not force_unique:
            obj.sort_time()

        return obj

    def _assert_range(self, valrange):
        assert valrange[0] <= valrange[1], \
            'Range must be in increasing order: (lo, hi)'
        return valrange
    
    def __repr__(self):
        s = '<EventList: {0} events;\n'.format(self.size)
        s += ' time range {0};\n channel range {1}>'.format(self.time_range,
                                                            self.channel_range)
        return s


class Bins():
    """A primitive class defining a set of histogram bins
    
    Parameters:
        counts (np.array): The array of counts in each bin
        lo_edges (np.array): The low-value edges of the bins
        hi_edges (np.array): The high-value edges of the bins
    """
    def __init__(self, counts, lo_edges, hi_edges):
        
        try:
            iter(counts)
            self._counts = np.asarray(counts).flatten()
        except:
            raise TypeError('counts must be an iterable')

        try:
            iter(lo_edges)
            self._lo_edges = np.asarray(lo_edges).flatten()
        except:
            raise TypeError('lo_edges must be an iterable')

        try:
            iter(hi_edges)
            self._hi_edges = np.asarray(hi_edges).flatten()
        except:
            raise TypeError('hi_edges must be an iterable')

        if (self._counts.size != self._lo_edges.size) or \
           (self._lo_edges.size != self._hi_edges.size):
            raise ValueError('counts, lo_edges, and hi_edges must all be ' \
                             'the same length')        

    @property
    def centroids(self):
        """(np.array): The centroids of the bins"""
        return (self.hi_edges + self.lo_edges) / 2.0

    @property
    def counts(self):
        """(np.array): The counts in each bin"""
        return self._counts

    @property
    def count_uncertainty(self):
        """(np.array): The count uncertainty in each bin"""
        return np.sqrt(self._counts)

    @property
    def hi_edges(self):
        """(np.array): The high-value edges of the bins"""
        return self._hi_edges

    @property
    def lo_edges(self):
        """(np.array): The low-value edges of the bins"""
        return self._lo_edges

    @property
    def range(self):
        """(float, float): The range of the bin edges"""
        if self.size > 0:
            return (self.lo_edges[0], self.hi_edges[-1])

    @property
    def rates(self):
        """(np.array): The count rate of each bin; counts/width"""
        return self.counts / self.widths

    @property
    def rate_uncertainty(self):
        """(np.array): The count rate uncertainty of each bin"""
        return self.count_uncertainty / self.widths

    @property
    def size(self):
        """(int): Number of bins"""
        return self.lo_edges.size

    @property
    def widths(self):
        """(np.array): The widths of the bins"""
        return self.hi_edges - self.lo_edges

    def closest_edge(self, val, which='either'):
        """Return the closest bin edge
        
        Args:
            val (float): Input value
            which (str, optional): Options are: 
                
                * 'either' - closest edge to val; 
                * 'low' - closest edge lower than val; or 
                * 'high' - closest edge higher than val. Default is 'either'
        
        Returns:        
            (float)
        """
        edges = np.concatenate((self.lo_edges, [self.hi_edges[-1]]))
        idx = np.argmin(np.abs(val - edges))
        if which == 'low' and (idx - 1) >= 0:
            if edges[idx] > val:
                idx -= 1
        elif (which == 'high') and (idx + 1) < edges.size:
            if edges[idx] < val:
                idx += 1
        else:
            pass
        return edges[idx]

    def slice(self, lo_edge, hi_edge):
        """Perform a slice over the range of the bins and return a new Bins 
        object. Note that lo_edge and hi_edge values that fall inside a bin will
        result in that bin being included.
        
        Args:
            lo_edge (float): The start of the slice
            hi_edge (float): The end of the slice
        
        Returns:
            (:class:`Bins`)
        """
        lo_snap = self.closest_edge(lo_edge, which='low')
        hi_snap = self.closest_edge(hi_edge, which='high')
        if lo_snap == hi_snap:
            mask = (self.lo_edges < hi_snap) & (self.hi_edges >= lo_snap)
        else:
            mask = (self.lo_edges < hi_snap) & (self.hi_edges > lo_snap)
        obj = Bins(self.counts[mask], self.lo_edges[mask], self.hi_edges[mask])
        return obj

    def __repr__(self):
        s = '<Bins: {0} bins;\n'.format(self.size)
        s += ' range {0}>'.format(self.range)
        return s
        

class ExposureBins(Bins):
    """A class defining a set of bins containing exposure information.

    Parameters:
        counts (np.array): The array of counts in each bin
        lo_edges (np.array): The low-value edges of the bins
        hi_edges (np.array): The high-value edges of the bins
        exposure (np.array): The exposure of each bin
        precalc_good_segments (bool, optional): If True, calculates contiguous
                                                bin segments on initialization.
                                                Default is True.
    """
    def __init__(self, counts, lo_edges, hi_edges, exposure, 
                 precalc_good_segments=True):
        super().__init__(counts, lo_edges, hi_edges)

        try:
            iter(exposure)
            self._exposure = np.asarray(exposure).flatten()
        except:
            raise TypeError('exposure must be an iterable')

        if (self._counts.size != self._exposure.size):
            raise ValueError('exposure must be the same size as counts')        
        
        self._good_segments = None
        if precalc_good_segments:
            self._good_segments = self._calculate_good_segments()

    @property
    def exposure(self):
        """(np.array): The exposure of each bin"""
        return self._exposure

    @property
    def rates(self):
        """(np.array): The count rate of each bin: counts/exposure"""
        r = np.zeros_like(self.exposure)
        mask = (self.exposure > 0.0)
        r[mask] = self.counts[mask] / self.exposure[mask]
        return r

    @property
    def rate_uncertainty(self):
        """(np.array): The count rate uncertainty of each bin"""
        r = np.zeros_like(self.exposure)
        mask = (self.exposure > 0.0)
        r[mask] = self.count_uncertainty[mask] / self.exposure[mask]
        return r
    
    def contiguous_bins(self):
        """Return a list of ExposureBins, each one containing a continuous 
        segment of data.  This is done by comparing the edges of each bin, and 
        if there is a gap between edges, the data is split into separate 
        ExposureBin objects, each containing a contiguous set of data.
        
        Returns
            (list of :class:`ExposureBins`)
        """
        if self._good_segments is not None:
            good_segments = self._good_segments
        else:
            good_segments = self._calculate_good_segments()
        bins = [self.slice(seg[0], seg[1]) for seg in good_segments]
        return bins

    def rebin(self, method, *args, tstart=None, tstop=None):
        """Rebin the ExposureBins object in given a binning function and return 
        a new ExposureBins object 
        
        The binning function should take as input an array of counts, 
        array of exposures, and an array of bin edges. Additional arguments 
        specific to the function are allowed. The function should return an 
        array of the new counts, new exposure, and new edges.
        
        Args:
            method (<function>): A binning function
            *args: Arguments to be passed to the binning function
            tstart (float, optional): 
                If set, defines the start time of the ExposureBins to be binned, 
                otherwise binning will begin at the time of the first bin edge.
            tstop (float, optional): 
                If set, defines the end time of the ExposureBins to be binned, 
                otherwise binning will end at the time of the last bin edge.        
        
        Returns:
            (:class:`ExposureBins`)
        """

        # empty bins
        empty = self.__class__([], [], [], [])

        # set the start and stop of the rebinning segment
        trange = self.range
        if tstart is None:
            tstart = trange[0]
        if tstop is None:
            tstop = trange[1]
        if tstart < trange[0]:
            tstart = trange[0]
        if tstop > trange[1]:
            tstop = trange[1]
        
        
        bins = self.contiguous_bins()
        new_histos = []
        for bin in bins:

            trange = bin.range

            # split the histogram into pieces so that we only rebin the piece
            # that needs to be rebinned
            pre = empty
            post = empty
            if (tstop < trange[0]) or (tstart > trange[1]):
                histo = empty
            elif tstop == trange[1]:
                if tstart > trange[0]:
                    pre = bin.slice(trange[0], self.closest_edge(tstart, 
                                                                 which='low'))
                histo = bin.slice(self.closest_edge(tstart, which='low'),
                                  trange[1])
            elif tstart == trange[0]:
                histo = bin.slice(trange[0],
                                  self.closest_edge(tstop, which='high'))
                if tstop < trange[1]:
                    post = bin.slice(self.closest_edge(tstop, which='high'), 
                                     trange[1])

            elif (tstart > trange[0]) or (tstop < trange[1]):
                pre = bin.slice(trange[0], 
                                self.closest_edge(tstart, which='low'))
                histo = bin.slice(self.closest_edge(tstart, which='low'),
                                  self.closest_edge(tstop, which='high'))
                post = bin.slice(self.closest_edge(tstop, which='high'), 
                                 trange[1])
            else:
                histo = bin

            # perform the rebinning and create a new histo with the 
            # rebinned rates
            if histo.size > 0:
                edges = np.append(histo.lo_edges, histo.hi_edges[-1])
                new_counts, new_exposure, new_edges = method(histo.counts,
                                                             histo.exposure,
                                                             edges, *args)
                new_histo = self.__class__(new_counts, new_edges[:-1],
                                           new_edges[1:],
                                           new_exposure)
            else:
                new_histo = bin

            # now merge the split histo back together again
            histos_to_merge = [i for i in (pre, new_histo, post) if
                               i.size > 0]
            
            if len(histos_to_merge) > 0: 
                new_histos.append(self.__class__.merge(histos_to_merge))

        new_histo = self.__class__.merge(new_histos)
        return new_histo
    
    def slice(self, tstart, tstop):
        """Perform a slice over a range and return a new ExposureBins object. 
        Note that tstart and tstop values that fall inside a bin will result in 
        that bin being included.

        Args:
            tstart (float): The start of the slice
            tstop (float): The end of the slice
        
        Returns:
            (:class:`ExposureBins`)
        """
        tstart_snap = self.closest_edge(tstart, which='low')
        tstop_snap = self.closest_edge(tstop, which='high')

        mask = (self.lo_edges < tstop_snap) & (self.hi_edges > tstart_snap)
        
        obj = self.__class__(self.counts[mask], self.lo_edges[mask],
                             self.hi_edges[mask], self.exposure[mask])
        return obj

    @classmethod
    def merge(cls, histos, **kwargs):
        """Merge multiple ExposureBins together.  Note that the ExposureBins
        to merged must be strictly non-overlapping.  The end of the range for
        one ExposureBins may be equal to the start of the range for another, or
        there may be a gap between the end of one or beginning of another,
        but there cannot be an overlap.  If you need to sum the counts between
        multiple ExposureBins with the same bin edges, see the 
        :meth:`~ExposureBins.sum` method.
        
        Args:
            histos (list of :class:`ExposureBins`): 
                A list containing the ExposureBins to be merged
        
        Returns:        
            (:class:`ExposureBins`)
        """
        num = len(histos)

        # sort by start time
        tstarts = np.concatenate([[histo.lo_edges[0]] for histo in histos])
        idx = np.argsort(tstarts)
        
        # concatenate the histos in order
        counts = histos[idx[0]].counts
        lo_edges = histos[idx[0]].lo_edges
        hi_edges = histos[idx[0]].hi_edges
        exposure = histos[idx[0]].exposure
        for i in range(1, num):
            bin_starts = histos[idx[i]].lo_edges
            # make sure there is no overlap
            mask = (bin_starts >= hi_edges[-1])
            if (~mask).sum() > 0:
                raise ValueError('Overlapping bins cannot be merged.  Only' \
                                 'non-overlapping bins can be merged.')
            
            counts = np.concatenate((counts, histos[idx[i]].counts[mask]))
            lo_edges = np.concatenate(
                (lo_edges, histos[idx[i]].lo_edges[mask]))
            hi_edges = np.concatenate(
                (hi_edges, histos[idx[i]].hi_edges[mask]))
            exposure = np.concatenate(
                (exposure, histos[idx[i]].exposure[mask]))

        # new ExposureBins object
        merged_bins = cls(counts, lo_edges, hi_edges, exposure, **kwargs)
        return merged_bins

    @classmethod
    def sum(cls, histos):
        """Sum multiple ExposureBins together if they have the same bin edges.
        If the exposures are different between the histograms, they will be 
        averaged.
        
        Args:
            histos (list of :class:`ExposureBins`):  
                A list containing the ExposureBins to be summed
        
        Returns:        
            (:class:`ExposureBins`)
        """
        counts = np.zeros(histos[0].size)
        for histo in histos:
            assert histo.size == histos[0].size, \
                "The histograms must all have the same size"
            assert np.all(histo.lo_edges == histos[0].lo_edges), \
                "The histograms must all have the same support"
            counts += histo.counts

        # averaged exposure
        exposure = np.mean([histo.exposure for histo in histos], axis=0)

        sum_bins = cls(counts, histos[0].lo_edges, histos[0].hi_edges,
                       exposure)
        return sum_bins

    def _calculate_good_segments(self):
        """Calculates the ranges of data that are contiguous segments
        
        Returns:
            ([(float, float), ...])
        """
        mask = (self.lo_edges[1:] != self.hi_edges[:-1])
        if mask.sum() == 0:
            return [self.range]
        times = np.concatenate(([self.lo_edges[0]], self.hi_edges[:-1][mask],
                                self.lo_edges[1:][mask], [self.hi_edges[-1]]))
        times.sort()
        return times.reshape(-1, 2).tolist()

    def __repr__(self):
        s = '<{0}: {1} bins;\n'.format(self.__class__.__name__, self.size)
        s += ' range {0};\n'.format(self.range)
        if self._good_segments is not None:
            s += ' {0} contiguous segments>'.format(len(self._good_segments))
        return s


class ChannelBins(ExposureBins):
    """A class defining a set of Energy Channel bins.
    """    
    @property
    def chan_nums(self):
        """(np.array): The channel numbers"""
        if self.size > 0:
            return self.lo_edges
            
    @classmethod
    def create(cls, counts, chan_nums, exposure, **kwargs):
        """Create a :class:`ChannelBins` object from a list of channel numbers.
        
        Args:
            counts (np.array): The array of counts in each bin
            chan_nums (np.array): The energy channel numbers
            exposure (np.array): The exposure of each bin
            precalc_good_segments (bool, optional): If True, calculates contiguous
                                                    bin segments on initialization.
                                                    Default is True.
        
        Returns:
            (:class:`ChannelBins`)
        """
        try:
            iter(exposure)
        except:
            counts = np.asarray(counts)
            exposure = np.full(counts.shape, exposure)
        
        chan_nums = np.asarray(chan_nums)
        return cls(counts, chan_nums, chan_nums + 1, exposure, **kwargs)
    
    @property
    def range(self):
        """(int, int): The channel range"""
        if self.size > 0:
            return (self.chan_nums[0], self.chan_nums[-1]) 

    def rebin(self, method, *args, chan_min=None, chan_max=None):        
        """Rebin the ChannelBins object given a binning function and return 
        a new ChannelBins object 
        
        The binning function should take as input an array of counts, 
        array of exposures, and an array of bin edges. Additional arguments 
        specific to the function are allowed. The function should return an 
        array of the new counts, new exposure, and new edges.
        
        Note::
          Edges for energy channels are treated as [chan_num, chan_num+1]
        
        Args:
            method (<function>): A binning function
            *args: Arguments to be passed to the binning function
            chan_min (float, optional): 
                If set, defines the minimum channel number of the ChannelBins to 
                be binned, otherwise binning will begin at the first channel 
                number.
            chan_max (float, optional): 
                If set, defines the maximum channel of the ChannelBins to be 
                binned, otherwise binning will end at the last channel number.        
        
        Returns:
            (:class:`ChannelBins`)
        """

        # empty bins
        empty = self.__class__([], [], [], [])

        # set the start and stop of the rebinning segment
        crange = self.range
        if chan_min is None:
            chan_min = crange[0]
        if chan_max is None:
            chan_max = crange[1]
        if chan_min < crange[0]:
            chan_min = crange[0]
        if chan_max > crange[1]:
            chan_max = crange[1]
        
        bins = self.contiguous_bins()
        new_histos = []
        for bin in bins:

            crange = bin.range
            # split the histogram into pieces so that we only rebin the piece
            # that needs to be rebinned
            pre = empty
            post = empty
            if (chan_max < crange[0]) or (chan_min > crange[1]):
                histo = empty
            
            elif chan_max == crange[1]:
                if chan_min > crange[0]:
                    pre = bin.slice(crange[0], chan_min-1)
                histo = bin.slice(chan_min, crange[1])
            
            elif chan_min == crange[0]:
                histo = bin.slice(crange[0], chan_max)
                if chan_max < crange[1]:
                    post = bin.slice(chan_max+1, crange[1])

            elif (chan_min > crange[0]) or (chan_max < crange[1]):
                pre = bin.slice(crange[0], chan_min-1)
                histo = bin.slice(chan_min, chan_max)
                post = bin.slice(chan_max+1, crange[1])
            
            else:
                histo = bin

            # perform the rebinning and create a new histo with the 
            # rebinned rates
            if histo.size > 0:
                edges = np.append(histo.lo_edges, histo.hi_edges[-1])
                new_counts, new_exposure, new_edges = method(histo.counts,
                                                             histo.exposure,
                                                             edges, *args)
                new_histo = self.__class__(new_counts, new_edges[:-1],
                                           new_edges[1:],
                                           new_exposure)
            else:
                new_histo = bin

            # now merge the split histo back together again
            histos_to_merge = [i for i in (pre, new_histo, post) if
                               i.size > 0]
            
            if len(histos_to_merge) > 0: 
                new_histos.append(self.__class__.merge(histos_to_merge))

        new_histo = self.__class__.merge(new_histos)
        return new_histo

    def slice(self, chan_min, chan_max):
        """Perform a slice over the range of the bins and return a new 
        ChannelBins object.
        
        Args:
            lo_edge (float): The start of the slice
            hi_edge (float): The end of the slice
        
        Returns:
            (:class:`ChannelBins`)
        """
        mask = (self.chan_nums <= chan_max) & (self.chan_nums >= chan_min)
        obj = self.__class__.create(self.counts[mask], self.chan_nums[mask], 
                                    self.exposure[mask])
        return obj


class TimeBins(ExposureBins):
    """A class defining a set of Time History (lightcurve) bins.

    Parameters:
        counts (np.array): The array of counts in each bin
        lo_edges (np.array): The low-value edges of the bins
        hi_edges (np.array): The high-value edges of the bins
        exposure (np.array): The exposure of each bin
    """
    def __init__(self, counts, lo_edges, hi_edges, exposure, **kwargs):
        super().__init__(counts, lo_edges, hi_edges, exposure, **kwargs)


class EnergyBins(ExposureBins):
    """A class defining a set of Energy (count spectra) bins.

    Parameters:
        counts (np.array): The array of counts in each bin
        lo_edges (np.array): The low-value edges of the bins
        hi_edges (np.array): The high-value edges of the bins
        exposure (np.array): The exposure of each bin
        precalc_good_segments (bool, optional): If True, calculates contiguous
                                                bin segments on initialization.
                                                Default is True.
    """
    def __init__(self, counts, lo_edges, hi_edges, exposure, **kwargs):
        try:
            iter(exposure)
        except:
            counts = np.asarray(counts)
            exposure = np.full(counts.shape, exposure)
        super().__init__(counts, lo_edges, hi_edges, exposure, **kwargs)

    @property
    def centroids(self):
        """(np.array): The centroids of the bins"""
        return np.sqrt(self.hi_edges * self.lo_edges)

    @property
    def rates_per_kev(self):
        """(np.array): Differential count rate"""
        return self.rates / self.widths
    
    @property
    def rate_uncertainty_per_kev(self):
        """(np.array): Differential rate uncertainty"""
        return self.rate_uncertainty / self.widths

    def rebin(self, method, *args, emin=None, emax=None):
        """Rebin the EnergyBins object in given a binning function and return a
        a new EnergyBins object 

        The binning function should take as input an array of counts, 
        array of exposures, and an array of bin edges. Additional arguments 
        specific to the function are allowed. The function should return an 
        array of the new counts, new exposure, and new edges.

        Args:
            method (<function>): A binning function
            *args: Arguments to be passed to the binning function
            emin (float, optional): 
                If set, defines the starting energy of the EnergyBins to be 
                binned, otherwise binning will begin at the first bin edge.
            emax (float, optional): 
                If set, defines the ending energy of the EnergyBins to be binned, 
                otherwise binning will end at the last bin edge.        
        
        Returns:
            (:class:`EnergyBins`)
        """
        histo = super().rebin(method, *args, tstart=emin, tstop=emax)
        histo._exposure = self.exposure[:histo.size]
        return histo

    @classmethod
    def sum(cls, histos):
        """Sum multiple EnergyBins together if they have the same energy range
        (support).  Example use would be summing two count spectra.
        
        Args:
            histos (list of :class:`EnergyBins`):  
                A list containing the EnergyBins to be summed
        
        Returns:        
            (:class:`EnergyBins`)
        """
        counts = np.zeros(histos[0].size)
        exposure = 0.0
        for histo in histos:
            assert histo.size == histos[0].size, \
                "The histograms must all have the same size"
            assert np.all(histo.lo_edges == histos[0].lo_edges), \
                "The histograms must all have the same support"
            counts += histo.counts
            exposure += histo.exposure

        sum_bins = cls(counts, histos[0].lo_edges, histos[0].hi_edges,
                       exposure)
        return sum_bins


class TimeChannelBins():
    """A class defining a set of 2D Time/Energy-Channel bins.

    Parameters:
        counts (np.array): The array of counts in each bin
        tstart (np.array): The low-value edges of the time bins
        tstop (np.array): The high-value edges of the time bins
        exposure (np.array): The exposure of each bin
        chan_nums (np.array): The channel numbers in ascending order
        quality (np.array, optional): The spectrum quality flag
        precalc_good_segments (bool, optional): If True, calculates the 
                                                good time and channel segments 
                                                on initialization.    
    """
    def __init__(self, counts, tstart, tstop, exposure, chan_nums,
                 quality=None, precalc_good_segments=True):

        try:
            iter(counts)
            self._counts = np.asarray(counts)
        except:
            raise TypeError('counts must be an iterable')
        if self._counts.ndim != 2:
            raise TypeError('counts must be a 2-dimensional array')
        
        try:
            iter(tstart)
            self._tstart = np.asarray(tstart).flatten()
        except:
            raise TypeError('tstart must be an iterable')
        try:
            iter(tstop)
            self._tstop = np.asarray(tstop).flatten()
        except:
            raise TypeError('tstop must be an iterable')
        try:
            iter(exposure)
            self._exposure = np.asarray(exposure).flatten()
        except:
            raise TypeError('exposure must be an iterable')

        try:
            iter(chan_nums)
            self._chan_nums = np.asarray(chan_nums, dtype=int).flatten()
        except:
            raise TypeError('chan_nums must be an iterable')
                
        if (self._tstart.size != self._tstop.size) or \
           (self._exposure.size != self._tstart.size):
            raise ValueError('tstart, tstop, and exposure must all be ' \
                             'the same length')        
        
        if (self._counts.shape[0] != self._tstart.size) or \
           (self._counts.shape[1] != self._chan_nums.size):
            raise ValueError('counts axis 0 must have same length as tstart, '\
                             'tstop, and exposure.  counts axis 1 must have '\
                             'same length as chan_nums.')
        
        if quality is not None:
            try:
                iter(quality)
                self._quality = np.asarray(quality).flatten()
            except:
                raise TypeError('quality must be an iterable')
            if self._quality.size != self._tstart.size:
                raise ValueError('quality must be same length as tstart')
        else:
            self._quality = np.zeros_like(self._tstart)
        
        
        self._good_time_segments = None
        self._good_channel_segments = None
        if (self.num_times > 0) and precalc_good_segments:
            self._good_time_segments = self._calculate_good_segments(
                                                                    self.tstart,
                                                                     self.tstop)
                
        if (self.num_chans > 0) and precalc_good_segments:
            chans0 = self.chan_nums
            chans1 = self.chan_nums + 1
            self._good_channel_segments = self._calculate_good_segments(chans0,
                                                                        chans1)

    @property
    def chan_nums(self):
        """(np.array): The channel numbers"""
        return self._chan_nums
    
    @property
    def channel_range(self):
        """(int, int): The channel number range"""
        if self.num_chans > 0:
            return (self.chan_nums[0], self.chan_nums[-1])
    
    @property
    def counts(self):
        """(np.array): The array of counts in each bin"""
        return self._counts

    @property
    def count_uncertainty(self):
        """ (np.array): The counts uncertainty in each bin"""
        return np.sqrt(self.counts)
    
    @property
    def exposure(self):
        """(np.array): The exposure of each bin"""
        return self._exposure
    
    @property
    def num_chans(self):
        """(int): The number of energy channels along the energy axis"""
        return self._chan_nums.size

    @property
    def num_times(self):
        """(int): The number of bins along the time axis"""
        return self._exposure.size

    @property
    def quality(self):
        """(np.array): The spectrum quality flag"""
        return self._quality

    @property
    def rates(self):
        """(np.array): The rates in each Time-Channel Bin"""
        return self.counts / (self.exposure[:, np.newaxis])

    @property
    def rate_uncertainty(self):
        """(np.array): The rate uncertainty in each bin"""
        return self.count_uncertainty / (self.exposure[:, np.newaxis])

    @property
    def size(self):
        """(int, int): The number of bins along both axes (num_times, num_chans)"""
        return self.counts.shape

    @property
    def time_centroids(self):
        """(np.array): The bin centroids along the time axis"""
        return (self.tstop + self.tstart) / 2.0

    @property
    def time_range(self):
        """(float, float): The range of the data along the time axis"""
        if self.num_times > 0:
            return (self.tstart[0], self.tstop[-1])

    @property
    def time_widths(self):
        """(np.array): The bin widths along the time axis"""
        return (self.tstop - self.tstart)

    @property
    def tstart(self):
        """(np.array): The low-value edges of the time bins"""
        return self._tstart
    
    @property
    def tstop(self):
        """(np.array): The high-value edges of the time bins"""
        return self._tstop

    def apply_ebounds(self, ebounds):
        """Apply an energy bounds calibration and return a TimeEnergyBins 
        object.
        
        Args:
            ebounds (:class:`Ebounds`): The energy bounds. Must match the number
                                        of channels as this object.
        
        Returns:
            (:class:`TimeEnergyBins`)
        """
        if not isinstance(ebounds, Ebounds):
            raise TypeError('ebounds must be an Ebounds object')
        
        if ebounds.num_intervals != self.num_chans:
            raise ValueError('ebounds must have the same number of channels ' \
                             'as this object.')
        
        return TimeEnergyBins(self.counts, self.tstart, self.tstop, 
                              self.exposure, ebounds.low_edges(),
                              ebounds.high_edges())
    
    def closest_time_edge(self, val, which='either'):
        """Return the closest time bin edge
        
        Args:
            val (float): Input value
            which (str, optional): Options are: 
                
                * 'either' - closest edge to val; 
                * 'low' - closest edge lower than val; 
                * 'high' - closest edge higher than val. Default is 'either'
        
        Returns:           
            (float)
        """
        edges = np.concatenate((self.tstart, [self.tstop[-1]]))
        return self._closest_edge(edges, val, which=which)

    def contiguous_channel_bins(self):
        """Return a list of TimeChannelBins, each one containing a contiguous
        channel segment of data.  This is done by comparing adjacent channel 
        numbers, and if there is a gap between channels, the data is split into 
        separate TimeChannelBins objects, each containing a 
        channel-contiguous set of data.
        
        Returns:
            (list of :class:`TimeChannelBins`)
        """
        if self._good_channel_segments is None:
            chans0 = self.chan_nums
            chans1 = self.chan_nums + 1
            good_segments = self._calculate_good_segments(chans0, chans1)
        else:
            good_segments = self._good_channel_segments
        bins = [self.slice_channels(seg[0], seg[1]) for seg in good_segments]
        return bins
    
    def contiguous_time_bins(self):
        """Return a list of TimeChannelBins, each one containing a contiguous
        time segment of data.  This is done by comparing the edges of each time
        bin, and if thereis a gap between edges, the data is split into 
        separate TimeChannelBins objects, each containing a time-contiguous set 
        of data.
        
        Returns:
            (list of :class:`TimeChannelBins`)
        """
        if self._good_time_segments is None:
            good_segments = self._calculate_good_segments(self.tstart,
                                                          self.tstop)
        else:
            good_segments = self._good_time_segments

        bins = [self.slice_time(seg[0], seg[1]) for seg in good_segments]
        return bins

    def get_exposure(self, time_ranges=None, scale=False):
        """Calculate the total exposure of a time range or time ranges of data
        
        Args:
            time_ranges ([(float, float), ...], optional): 
                The time range or time ranges over which to calculate the 
                exposure. If omitted, calculates the total exposure of the data        
            scale (bool, optional): 
                If True and the time ranges don't match up with the data binning, 
                will scale the exposure based on the requested time range. 
                Default is False.
        
        Returns:
            (float)
        """
        if time_ranges is None:
            time_ranges = [self.time_range]
        try:
            iter(time_ranges[0])
        except:
            time_ranges = [time_ranges]
        exposure = 0.0

        for i in range(len(time_ranges)):
            mask = self._slice_time_mask(*self._assert_range(time_ranges[i]))
            dt = (time_ranges[i][1] - time_ranges[i][0])
            data_exp = np.sum(self.exposure[mask])
            dts = np.sum(self.tstop[mask] - self.tstart[mask])
            if dts > 0:
                if scale:
                    exposure += data_exp * (dt / dts)
                else:
                    exposure += data_exp

        return exposure

    def integrate_channels(self, chan_min=None, chan_max=None):
        """Integrate the histogram over the channel axis (producing a 
        lightcurve). Limits on the integration smaller than the full range can 
        be set.
        
        Args:
            chan_min (int, optional): 
                The low end of the integration range. If not set, uses the 
                lowest energy channel of the histogram
            chan_max (int, optional): 
                The high end of the integration range. If not set, uses the 
                highest energy channel of the histogram
        
        Returns:           
            (:class:`TimeBins`)
        """
        if chan_min is None:
            chan_min = self.channel_range[0]
        if chan_max is None:
            chan_max = self.channel_range[1]

        mask = (self.chan_nums <= chan_max) & (self.chan_nums >= chan_min)
        counts = self.counts[:, mask].sum(axis=1)

        obj = TimeBins(counts, self.tstart, self.tstop, self.exposure)
        return obj

    def integrate_time(self, tstart=None, tstop=None):
        """Integrate the histogram over the time axis (producing an energy 
        channel spectrum). Limits on the integration smaller than the full range 
        can be set.
        
        Args:
            tstart (float, optional): 
                The low end of the integration range. If not set, uses the 
                lowest time edge of the histogram
            tstop (float, optional): 
                The high end of the integration range. If not set, uses the 
                highest time edge of the histogram
        
        Returns:           
            (:class:`ChannelBins`)
        """
        if tstart is None:
            tstart = self.time_range[0]
        if tstop is None:
            tstop = self.time_range[1]

        mask = self._slice_time_mask(tstart, tstop)
        counts = np.sum(self.counts[mask, :], axis=0)
        exposure = np.sum(self.exposure[mask])
        exposure = np.full(counts.size, exposure)

        obj = ChannelBins.create(counts, self.chan_nums, exposure)
        return obj

    def rebin_channels(self, method, *args, chan_min=None, chan_max=None):
        """Rebin the TimeChannelBins object along the energy axis given a 
        binning function and return a new TimeChannelBins object.
        
        The binning function should take as input an array of counts, 
        array of exposures, and an array of bin edges. Additional arguments 
        specific to the function are allowed. The function should return an 
        array of the new counts, new exposure, and new edges.
        
        Args:
            method (<function>):  A binning function
            *args:  Arguments to be passed to the binning function
            chan_min (int, optional): 
                If set, defines the starting channel of the TimeChannelBins 
                to be binned, otherwise binning will begin at the the first 
                channel.
            chan_max (int, optional): 
                If set, defines the ending channel of the TimeChannelBins to 
                be binned, otherwise binning will end at the last channel.
        
        Returns:        
            (:class:`TimeChannelBins`)
        """
        # empty bins
        empty = self.__class__(np.array([[]]).reshape(0,0), [], [], [], [])

        # set the start and stop of the rebinning segment
        chan_range = self.channel_range
        if chan_min is None:
            chan_min = chan_range[0]
        if chan_max is None:
            chan_max = chan_range[1]
        if chan_min < chan_range[0]:
            chan_min = chan_range[0]
        if chan_max > chan_range[1]:
            chan_max = chan_range[1]

        bins = self.contiguous_channel_bins()
        new_histos = []
        for bin in bins:
            # split the histogram into pieces so that we only rebin the piece
            # that needs to be rebinned
            pre = empty
            post = empty
            crange = bin.channel_range
            if (chan_max < crange[0]) or (chan_min > crange[1]):
                histo = empty
            
            elif chan_max == crange[1]:
                if chan_min > crange[0]:
                    pre = bin.slice_channels(crange[0], chan_min-1)
                histo = bin.slice_channels(chan_min, crange[1])
                
            elif chan_min == crange[0]:
                histo = bin.slice_channels(crange[0], chan_max)
                if chan_max < crange[1]:
                    post = bin.slice_channels(chan_max+1, crange[1])
            
            elif (chan_min > crange[0]) and (chan_max < crange[1]):
                pre = bin.slice_channels(crange[0], chan_min-1)
                histo = bin.slice_channels(chan_min, chan_max)
                post = bin.slice_channels(chan_max+1, crange[1])
            
            else:
                histo = bin
                
            # perform the rebinning and create a new histo with the 
            # rebinned rates
            if histo.num_chans > 0:
                edges = np.append(histo.chan_nums, histo.chan_nums[-1]+1)
                num_times, num_chans = histo.size
                new_counts = []
                for i in range(num_times):
                    exposure = np.full(num_chans, histo.exposure[i])
                    new_cts, _, new_edges = method(histo.counts[i, :],
                                                   exposure,
                                                   edges, *args)
                    new_counts.append(new_cts)
                new_counts = np.array(new_counts)
                new_histo = TimeChannelBins(new_counts, bin.tstart, 
                                               bin.tstop, bin.exposure, 
                                               new_edges[:-1])
            
            else:
                new_histo = bin
            
            # now merge the split histo back together again
            histos_to_merge = [i for i in (pre, new_histo, post) if
                               i.num_chans > 0]
            
            if len(histos_to_merge) > 0: 
                new_histos.append(TimeChannelBins.merge_channels(histos_to_merge))

        new_histo = TimeChannelBins.merge_channels(new_histos)

        return new_histo

    def rebin_time(self, method, *args, tstart=None, tstop=None):
        """Rebin the TimeChannelBins object along the time axis given a binning 
        function and return a new TimeChannelBins object.
        
        The binning function should take as input an array of counts, 
        array of exposures, and an array of bin edges. Additional arguments 
        specific to the function are allowed. The function should return an 
        array of the new counts, new exposure, and new edges.
        
        Args:
            method (<function>): A binning function
            *args:  Arguments to be passed to the binning function
            tstart (float, optional): 
                If set, defines the start time of the TimeChannelBins to be 
                binned, otherwise binning will begin at the time of the first 
                bin edge.
            tstop (float, optional): 
                If set, defines the end time of the TimeChannelBins to be 
                binned, otherwise binning will end at the time of the last 
                bin edge.
        
        Returns:       
            (:class:`TimeChannelBins`)
        """

        # empty bins
        empty = self.__class__(np.array([[]]).reshape(0,0), [], [], [], [], [])

        # set the start and stop of the rebinning segment
        trange = self.time_range
        if tstart is None:
            tstart = trange[0]
        if tstop is None:
            tstop = trange[1]
        if tstart < trange[0]:
            tstart = trange[0]
        if tstop > trange[1]:
            tstop = trange[1]

        bins = self.contiguous_time_bins()
        new_histos = []
        for bin in bins:
            trange = bin.time_range
            # split the histogram into pieces so that we only rebin the piece
            # that needs to be rebinned
            pre = empty
            post = empty
            if (tstop < trange[0]) or (tstart > trange[1]):
                histo = empty
            elif tstop == trange[1]:
                if tstart > trange[0]:
                    pre = bin.slice_time(trange[0], 
                                         self.closest_time_edge(tstart, 
                                                                which='low'))
                histo = bin.slice_time(self.closest_time_edge(tstart,
                                                              which='low'),
                                       trange[1])
            elif tstart == trange[0]:
                histo = bin.slice_time(trange[0],
                                       self.closest_time_edge(tstop,
                                                              which='high'))
                if tstop < trange[1]:
                    post = bin.slice_time(self.closest_time_edge(tstop, 
                                                                  which='high'), 
                                          trange[1])
                           
            elif (tstart > trange[0]) and (tstop < trange[1]):
                pre = bin.slice_time(trange[0], 
                                    self.closest_time_edge(tstart, which='low'))
                histo = bin.slice_time(
                                   self.closest_time_edge(tstart, which='low'),
                                   self.closest_time_edge(tstop, which='high'))
                post = bin.slice_time(self.closest_time_edge(tstop, which='high'), 
                                      trange[1])
            else:
                histo = bin
            
            # perform the rebinning and create a new histo with the 
            # rebinned rates
            if histo.num_times > 0:
                edges = np.append(histo.tstart, histo.tstop[-1])
                new_counts = []
                for i in range(bin.num_chans):
                    new_cts, new_exposure, new_edges = method(
                        histo.counts[:, i],
                        histo.exposure,
                        edges, *args)
                    new_counts.append(new_cts)
                new_counts = np.array(new_counts).T

                new_histo = TimeChannelBins(new_counts, new_edges[:-1],
                                           new_edges[1:], new_exposure, 
                                           bin.chan_nums)
            else:
                new_histo = bin

            # now merge the split histo back together again
            histos_to_merge = [i for i in (pre, new_histo, post) if
                               i.num_times > 0]
            
            if len(histos_to_merge) > 0: 
                new_histos.append(TimeChannelBins.merge_time(histos_to_merge))

        new_histo = TimeChannelBins.merge_time(new_histos)

        return new_histo

    def slice_channels(self, chan_min, chan_max):
        """Perform a slice over an energy range and return a new TimeChannelBins 
        object. Note that chan_min and chan_max values that fall inside a bin will 
        result in that bin being included.
        
        Args:
            emin (float): The start of the slice
            emax (float): The end of the slice
        
        Returns:           
            (:class:`TimeChannelBins`)
        """
        mask = (self.chan_nums <= chan_max) & (self.chan_nums >= chan_min)
        obj = TimeChannelBins(self.counts[:, mask], self.tstart, self.tstop,
                             self.exposure, self.chan_nums[mask],
                             quality=self.quality)
        return obj

    def slice_time(self, tstart, tstop):
        """Perform a slice over a time range and return a new TimeChannelBins 
        object. Note that tstart and tstop values that fall inside a bin will 
        result in that bin being included.
        
        Args:
            tstart (float): The start of the slice
            tstop (float): The end of the slice
        
        Returns:           
            (:class:`TimeChannelBins`)
        """
        mask = self._slice_time_mask(tstart, tstop)
        cls = type(self)
        obj = cls(self.counts[mask, :], self.tstart[mask], self.tstop[mask], 
                  self.exposure[mask], self.chan_nums, 
                  quality=self.quality[mask])
        return obj

    @classmethod
    def merge_channels(cls, histos, **kwargs):
        """Merge multiple TimeChannelBins together along the channel axis.
        
        Args:
            histos (list of :class:`TimeChannelBins`): 
                A list containing the TimeChannelBins to be merged
        
        Returns:
            (:class:`TimeChannelBins`)
        """
        num = len(histos)
        # sort by channel edge
        chan_mins = np.concatenate([[histo.chan_nums[0]] for histo in histos])
        idx = np.argsort(chan_mins)

        # concatenate the histos in order
        counts = histos[idx[0]].counts
        tstart = histos[idx[0]].tstart
        tstop = histos[idx[0]].tstop
        exposure = histos[idx[0]].exposure
        quality = histos[idx[0]].quality
        chan_nums = histos[idx[0]].chan_nums
        for i in range(1, num):
            chan_starts = histos[idx[i]].chan_nums
            # make sure there is no overlap
            mask = (chan_starts >= chan_nums[-1])
            if (~mask).sum() > 0:
                raise ValueError('Overlapping bins cannot be merged.  Only' \
                                 'non-overlapping bins can be merged.')

            counts = np.hstack((counts, histos[idx[i]].counts[:, mask]))
            chan_nums = np.concatenate((chan_nums, histos[idx[i]].chan_nums[mask]))

        # new TimeChannelBins object
        merged_bins = cls(counts, tstart, tstop, exposure, chan_nums, 
                          quality=quality, **kwargs)
        return merged_bins

    @classmethod
    def merge_time(cls, histos, **kwargs):
        """Merge multiple TimeChannelBins together along the time axis.
        
        Args:
            histos (list of :class:`TimeChannelBins`): 
                A list containing the TimeChannelBins to be merged
        
        Returns:
            (:class:`TimeChannelBins`)
        """
        num = len(histos)
        # sort by start time
        tstarts = np.concatenate([[histo.tstart[0]] for histo in histos])
        idx = np.argsort(tstarts)

        # concatenate the histos in order
        counts = histos[idx[0]].counts
        tstart = histos[idx[0]].tstart
        tstop = histos[idx[0]].tstop
        exposure = histos[idx[0]].exposure
        chan_nums = histos[idx[0]].chan_nums
        quality = histos[idx[0]].quality
        for i in range(1, num):
            bin_starts = histos[idx[i]].tstart
            # make sure there is no overlap
            mask = (bin_starts >= tstop[-1])
            if (~mask).sum() > 0:
                raise ValueError('Overlapping bins cannot be merged.  Only' \
                                 'non-overlapping bins can be merged.')

            counts = np.vstack((counts, histos[idx[i]].counts[mask, :]))
            tstart = np.concatenate((tstart, histos[idx[i]].tstart[mask]))
            tstop = np.concatenate((tstop, histos[idx[i]].tstop[mask]))
            exposure = np.concatenate(
                (exposure, histos[idx[i]].exposure[mask]))
            quality = np.concatenate((quality, histos[idx[i]].quality[mask]))

        # new TimeChannelBins object
        merged_bins = cls(counts, tstart, tstop, exposure, chan_nums, 
                          quality=quality, **kwargs)
        return merged_bins

    def _assert_range(self, valrange):
        assert valrange[0] <= valrange[1], \
            'Range must be in increasing order: (lo, hi)'
        return valrange
    
    def _calculate_good_segments(self, lo_edges, hi_edges):
        """Calculates the ranges of data that are contiguous segments
        
        Args:
            lo_edges (np.array): The lower bin edges
            hi_edges (np.array): The upper bin edges
        
        Returns:           
            ([(float, float), ...])
        """
        mask = (lo_edges[1:] != hi_edges[:-1])
        if mask.sum() == 0:
            return [(lo_edges[0], hi_edges[-1])]
        edges = np.concatenate(([lo_edges[0]], hi_edges[:-1][mask],
                                lo_edges[1:][mask], [hi_edges[-1]]))
        edges.sort()
        return edges.reshape(-1, 2).tolist()

    def _closest_edge(self, edges, val, which='either'):
        """Return the closest time bin edge
        
        Args:
            val (float): Input value
            which (str, optional): Options are: 
                
                * 'either' - closest edge to val; 
                * 'low' - closest edge lower than val; 
                * 'high' - closest edge higher than val. Default is 'either'
        
        Returns:           
            (float)
        """
        idx = np.argmin(np.abs(val - edges))
        if which == 'low':
            if (edges[idx] > val) and (idx - 1) >= 0:
                idx -= 1
        elif (which == 'high') and (idx + 1) < edges.size:
            if edges[idx] < val:
                idx += 1
        else:
            pass
        return edges[idx]

    def _slice_time_mask(self, tstart, tstop):
        tstart_snap = self.closest_time_edge(tstart, which='low')
        tstop_snap = self.closest_time_edge(tstop, which='high')
        mask = (self.tstart < tstop_snap) & (self.tstop > tstart_snap)
        return mask

    def __repr__(self):
        s = '<{0}: {1} time bins;\n'.format(self.__class__.__name__, 
                                            self.num_times)
        s += ' time range {0};\n'.format(self.time_range)
        if self._good_time_segments is not None:
            s += ' {0} time segments;\n'.format(len(self._good_time_segments))
        
        s += ' {0} channels;\n'.format(self.num_chans)
        s += ' channel range {0}'.format(self.channel_range)
        if self._good_channel_segments is not None:
            s += ';\n {0} channel segments'.format(len(self._good_channel_segments))
        
        return s+'>'


class TimeEnergyBins():
    """A class defining a set of 2D Time-Energy bins.

    Parameters:
        counts (np.array): The array of counts in each bin
        tstart (np.array): The low-value edges of the time bins
        tstop (np.array): The high-value edges of the time bins
        exposure (np.array): The exposure of each bin
        emin (np.array): The low-value edges of the energy bins
        emax (np.array): The high-value edges of the energy bins
        quality (np.array, optional): The spectrum quality flag
        precalc_good_segments (bool, optional): If True, calculates the 
                                                good time and energy segments on
                                                initialization.    
    """
    def __init__(self, counts, tstart, tstop, exposure, emin, emax,
                 quality=None, precalc_good_segments=True):

        try:
            iter(counts)
            self._counts = np.asarray(counts)
        except:
            raise TypeError('counts must be an iterable')
        if self._counts.ndim != 2:
            raise TypeError('counts must be a 2-dimensional array')
        
        try:
            iter(tstart)
            self._tstart = np.asarray(tstart).flatten()
        except:
            raise TypeError('tstart must be an iterable')
        try:
            iter(tstop)
            self._tstop = np.asarray(tstop).flatten()
        except:
            raise TypeError('tstop must be an iterable')
        try:
            iter(exposure)
            self._exposure = np.asarray(exposure).flatten()
        except:
            raise TypeError('exposure must be an iterable')

        try:
            iter(emin)
            self._emin = np.asarray(emin).flatten()
        except:
            raise TypeError('emin must be an iterable')
        try:
            iter(emax)
            self._emax = np.asarray(emax).flatten()
        except:
            raise TypeError('emax must be an iterable')
        
        if (self._tstart.size != self._tstop.size) or \
           (self._exposure.size != self._tstart.size):
            raise ValueError('tstart, tstop, and exposure must all be ' \
                             'the same length')        

        if (self._emin.size != self._emax.size):
            raise ValueError('emin and emax must be the same length')        
        
        if (self._counts.shape[0] != self._tstart.size) or \
           (self._counts.shape[1] != self._emin.size):
            raise ValueError('counts axis 0 must have same length as tstart, '\
                             'tstop, and exposure.  counts axis 1 must have '\
                             'same length as emin and emax.')
        
        if quality is not None:
            try:
                iter(quality)
                self._quality = np.asarray(quality).flatten()
            except:
                raise TypeError('quality must be an iterable')
            if self._quality.size != self._tstart.size:
                raise ValueError('quality must be same length as tstart')
        else:
            self._quality = np.zeros_like(self._tstart)
        
        
        self._good_time_segments = None
        self._good_energy_segments = None
        if (self.num_times > 0) and precalc_good_segments:
            self._good_time_segments = self._calculate_good_segments(
                                                                    self.tstart,
                                                                     self.tstop)
                
        if (self.num_chans > 0) and precalc_good_segments:
            self._good_energy_segments = self._calculate_good_segments(
                                                                      self.emin,
                                                                      self.emax)

    @property
    def chan_widths(self):
        """(np.array): The bin widths along the energy axis"""
        return (self.emax - self.emin)

    @property
    def counts(self):
        """(np.array): The array of counts in each bin"""
        return self._counts

    @property
    def count_uncertainty(self):
        """ (np.array): The counts uncertainty in each bin"""
        return np.sqrt(self.counts)
    
    @property
    def emax(self):
        """(np.array): The high-value edges of the energy bins"""
        return self._emax

    @property
    def emin(self):
        """(np.array): The low-value edges of the energy bins"""
        return self._emin
    
    @property
    def energy_centroids(self):
        """(np.array): The bin centroids along the energy axis"""
        return np.sqrt(self.emin * self.emax)

    @property
    def energy_range(self):
        """(float, float): The range of the data along the energy axis"""
        if self.num_chans > 0:
            return (self.emin[0], self.emax[-1])

    @property
    def exposure(self):
        """(np.array): The exposure of each bin"""
        return self._exposure
    
    @property
    def num_chans(self):
        """(int): The number of energy channels along the energy axis"""
        return self._emin.size

    @property
    def num_times(self):
        """(int): The number of bins along the time axis"""
        return self._exposure.size

    @property
    def quality(self):
        """(np.array): The spectrum quality flag"""
        return self._quality

    @property
    def rates(self):
        """(np.array): The rates in each Time-Energy Bin"""
        return self.counts / (self.exposure[:, np.newaxis])

    @property
    def rates_per_kev(self):
        """(np.array): The differential rates in units of counts/s/keV"""
        return self.rates / self.chan_widths[np.newaxis,:]

    @property
    def rate_uncertainty(self):
        """(np.array): The rate uncertainty in each bin"""
        return self.count_uncertainty / (self.exposure[:, np.newaxis])

    @property
    def rate_uncertainty_per_kev(self):
        """The differential rate uncertainty in units of counts/s/keV"""
        return self.rate_uncertainty / self.chan_widths[np.newaxis, :]

    @property
    def size(self):
        """(int, int): The number of bins along both axes (numtimes, num_chans)"""
        return self.counts.shape

    @property
    def time_centroids(self):
        """(np.array): The bin centroids along the time axis"""
        return (self.tstop + self.tstart) / 2.0

    @property
    def time_range(self):
        """(float, float): The range of the data along the time axis"""
        if self.num_times > 0:
            return (self.tstart[0], self.tstop[-1])

    @property
    def time_widths(self):
        """(np.array): The bin widths along the time axis"""
        return (self.tstop - self.tstart)

    @property
    def tstart(self):
        """(np.array): The low-value edges of the time bins"""
        return self._tstart
    
    @property
    def tstop(self):
        """(np.array): The high-value edges of the time bins"""
        return self._tstop

    def closest_energy_edge(self, val, which='either'):
        """Return the closest energy bin edge

        Args:
            val (float): Input value
            which (str, optional): Options are: 
                
                * 'either' - closest edge to val; 
                * 'low' - closest edge lower than val; 
                * 'high' - closest edge higher than val. Default is 'either'
        
        Returns:           
            (float)
        """
        edges = np.concatenate((self.emin, [self.emax[-1]]))
        return self._closest_edge(edges, val, which=which)

    def closest_time_edge(self, val, which='either'):
        """Return the closest time bin edge
        
        Args:
            val (float): Input value
            which (str, optional): Options are: 
                
                * 'either' - closest edge to val; 
                * 'low' - closest edge lower than val; 
                * 'high' - closest edge higher than val. Default is 'either'
        
        Returns:           
            (float)
        """
        edges = np.concatenate((self.tstart, [self.tstop[-1]]))
        return self._closest_edge(edges, val, which=which)

    def contiguous_energy_bins(self):
        """Return a list of TimeEnergyBins, each one containing a contiguous
        energy segment of data.  This is done by comparing the edges of each
        energy bin, and if thereis a gap between edges, the data is split into 
        separate TimeEnergyBin objects, each containing an energy-contiguous set 
        of data.
        
        Returns:
            (list of :class:`TimeEnergyBins`)
        """
        if self._good_energy_segments is None:
            good_segments = self._calculate_good_segments(self.emin, self.emax)
        else:
            good_segments = self._good_energy_segments
        bins = [self.slice_energy(seg[0], seg[1]) for seg in
                self._good_energy_segments]
        return bins

    def contiguous_time_bins(self):
        """Return a list of TimeEnergyBins, each one containing a contiguous
        time segment of data.  This is done by comparing the edges of each time
        bin, and if thereis a gap between edges, the data is split into 
        separate TimeEnergyBin objects, each containing a time-contiguous set 
        of data.
        
        Returns:
            (list of :class:`TimeEnergyBins`)
        """
        if self._good_time_segments is None:
            good_segments = self._calculate_good_segments(self.tstart,
                                                          self.tstop)
        else:
            good_segments = self._good_time_segments

        bins = [self.slice_time(seg[0], seg[1]) for seg in good_segments]
        return bins

    def get_exposure(self, time_ranges=None, scale=False):
        """Calculate the total exposure of a time range or time ranges of data
        
        Args:
            time_ranges ([(float, float), ...], optional): 
                The time range or time ranges over which to calculate the 
                exposure. If omitted, calculates the total exposure of the data        
            scale (bool, optional): 
                If True and the time ranges don't match up with the data binning, 
                will scale the exposure based on the requested time range. 
                Default is False.
        
        Returns:
            (float)
        """
        if time_ranges is None:
            time_ranges = [self.time_range]
        try:
            iter(time_ranges[0])
        except:
            time_ranges = [time_ranges]
        exposure = 0.0

        for i in range(len(time_ranges)):
            mask = self._slice_time_mask(*self._assert_range(time_ranges[i]))
            dt = (time_ranges[i][1] - time_ranges[i][0])
            data_exp = np.sum(self.exposure[mask])
            dts = np.sum(self.tstop[mask] - self.tstart[mask])
            if dts > 0:
                if scale:
                    exposure += data_exp * (dt / dts)
                else:
                    exposure += data_exp

        return exposure

    def integrate_energy(self, emin=None, emax=None):
        """Integrate the histogram over the energy axis (producing a lightcurve).
        Limits on the integration smaller than the full range can be set.
        
        Args:
            emin (float, optional): 
                The low end of the integration range. If not set, uses the 
                lowest energy edge of the histogram
            emax (float, optional): 
                The high end of the integration range. If not set, uses the 
                highest energy edge of the histogram
        
        Returns:           
            (:class:`TimeBins`)
        """
        if emin is None:
            emin = self.energy_range[0]
        if emax is None:
            emax = self.energy_range[1]

        mask = self._slice_energy_mask(emin, emax)
        counts = np.sum(self.counts[:, mask], axis=1)

        obj = TimeBins(counts, self.tstart, self.tstop, self.exposure)
        return obj

    def integrate_time(self, tstart=None, tstop=None):
        """Integrate the histogram over the time axis (producing a count rate
        spectrum). Limits on the integration smaller than the full range can 
        be set.
        
        Args:
            tstart (float, optional): 
                The low end of the integration range. If not set, uses the 
                lowest time edge of the histogram
            tstop (float, optional): 
                The high end of the integration range. If not set, uses the 
                highest time edge of the histogram
        
        Returns:           
            (:class:`EnergyBins`)
        """
        if tstart is None:
            tstart = self.time_range[0]
        if tstop is None:
            tstop = self.time_range[1]

        mask = self._slice_time_mask(tstart, tstop)
        counts = np.sum(self.counts[mask, :], axis=0)
        exposure = np.sum(self.exposure[mask])
        exposure = np.full(counts.size, exposure)

        obj = EnergyBins(counts, self.emin, self.emax, exposure)
        return obj

    def rebin_energy(self, method, *args, emin=None, emax=None):
        """Rebin the TimeEnergyBins object along the energy axis given a binning 
        function and return a new TimeEnergyBins object.
        
        The binning function should take as input an array of counts, 
        array of exposures, and an array of bin edges. Additional arguments 
        specific to the function are allowed. The function should return an 
        array of the new counts, new exposure, and new edges.
        
        Args:
            method (<function>):  A binning function
            *args:  Arguments to be passed to the binning function
            emin (float, optional): 
                If set, defines the starting edge of the TimeEnergyBins to be 
                binned, otherwise binning will begin at the the first bin edge.
            emax (float, optional): 
                If set, defines the ending edge of the TimeEnergyBins to be 
                binned, otherwise binning will end at the last bin edge.
        
        Returns:        
            (:class:`TimeEnergyBins`)
        """

        # empty bins
        empty = self.__class__(np.array([[]]).reshape(0,0), [], [], [], [], [])

        # set the start and stop of the rebinning segment
        erange = self.energy_range
        if emin is None:
            emin = erange[0]
        if emax is None:
            emax = erange[1]
        if emin < erange[0]:
            emin = erange[0]
        if emax > erange[1]:
            emax = erange[1]

        bins = self.contiguous_energy_bins()
        new_histos = []
        for bin in bins:
            # split the histogram into pieces so that we only rebin the piece
            # that needs to be rebinned
            pre = empty
            post = empty
            erange = bin.energy_range
            if (emax < erange[0]) or (emin > erange[1]):
                histo = empty
            elif emax == erange[1]:
                if emin > erange[0]:
                    pre = bin.slice_energy(erange[0], 
                                           self.closest_energy_edge(emin, 
                                                                   which='low'))
                histo = bin.slice_energy(self.closest_energy_edge(emin,
                                                                  which='low'),
                                         erange[1])
            elif emin == erange[0]:
                histo = bin.slice_energy(erange[0],
                                         self.closest_energy_edge(emax,
                                                                  which='high'))
                if emax < erange[1]:
                    post = bin.slice_energy(self.closest_energy_edge(emax, 
                                                                  which='high'), 
                                            erange[1])
            elif (emin > erange[0]) and (emax < erange[1]):
                pre = bin.slice_energy(erange[0], 
                                    self.closest_energy_edge(emin, which='low'))
                histo = bin.slice_energy(
                                  self.closest_energy_edge(emin, which='low'),
                                  self.closest_energy_edge(emax, which='high'))
                post = bin.slice_energy(self.closest_energy_edge(emax, 
                                        which='high'), erange[1])
            else:
                histo = bin
                
            # perform the rebinning and create a new histo with the 
            # rebinned rates
            if histo.num_chans > 0:
                edges = np.append(histo.emin, histo.emax[-1])
                num_times, num_chans = histo.size
                new_counts = []
                for i in range(num_times):
                    exposure = np.full(num_chans, histo.exposure[i])
                    new_cts, _, new_edges = method(histo.counts[i, :],
                                                   exposure,
                                                   edges, *args)
                    new_counts.append(new_cts)
                new_counts = np.array(new_counts)
                new_histo = TimeEnergyBins(new_counts, bin.tstart, bin.tstop,
                                           bin.exposure, new_edges[:-1],
                                           new_edges[1:])
            else:
                new_histo = bin
            
            # now merge the split histo back together again
            histos_to_merge = [i for i in (pre, new_histo, post) if
                               i.num_chans > 0]
            
            if len(histos_to_merge) > 0: 
                new_histos.append(TimeEnergyBins.merge_energy(histos_to_merge))

        new_histo = TimeEnergyBins.merge_energy(new_histos)

        return new_histo

    def rebin_time(self, method, *args, tstart=None, tstop=None):
        """Rebin the TimeEnergyBins object along the time axis given a binning 
        function and return a new TimeEnergyBins object 
        
        The binning function should take as input an array of counts, 
        array of exposures, and an array of bin edges. Additional arguments 
        specific to the function are allowed. The function should return an 
        array of the new counts, new exposure, and new edges.
        
        Args:
            method (<function>): A binning function
            *args:  Arguments to be passed to the binning function
            tstart (float, optional): 
                If set, defines the start time of the TimeEnergyBins to be 
                binned, otherwise binning will begin at the time of the first 
                bin edge.
            tstop (float, optional): 
                If set, defines the end time of the TimeEnergyBins to be 
                binned, otherwise binning will end at the time of the last 
                bin edge.
        
        Returns:       
            (:class:`TimeEnergyBins`)
        """

        # empty bins
        empty = self.__class__(np.array([[]]).reshape(0,0), [], [], [], [], [])

        # set the start and stop of the rebinning segment
        trange = self.time_range
        if tstart is None:
            tstart = trange[0]
        if tstop is None:
            tstop = trange[1]
        if tstart < trange[0]:
            tstart = trange[0]
        if tstop > trange[1]:
            tstop = trange[1]

        bins = self.contiguous_time_bins()
        new_histos = []
        for bin in bins:
            trange = bin.time_range
            # split the histogram into pieces so that we only rebin the piece
            # that needs to be rebinned
            pre = empty
            post = empty
            if (tstop < trange[0]) or (tstart > trange[1]):
                histo = empty
            elif tstop == trange[1]:
                if tstart > trange[0]:
                    pre = bin.slice_time(trange[0], 
                                         self.closest_time_edge(tstart, 
                                                                which='low'))
                histo = bin.slice_time(self.closest_time_edge(tstart,
                                                              which='low'),
                                       trange[1])
            elif tstart == trange[0]:
                histo = bin.slice_time(trange[0],
                                       self.closest_time_edge(tstop,
                                                              which='high'))
                if tstop < trange[1]:
                    post = bin.slice_time(self.closest_time_edge(tstop, 
                                                                  which='high'), 
                                          trange[1])
                           
            elif (tstart > trange[0]) and (tstop < trange[1]):
                pre = bin.slice_time(trange[0], 
                                    self.closest_time_edge(tstart, which='low'))
                histo = bin.slice_time(
                                   self.closest_time_edge(tstart, which='low'),
                                   self.closest_time_edge(tstop, which='high'))
                post = bin.slice_time(self.closest_time_edge(tstop, which='high'), 
                                      trange[1])
            else:
                histo = bin
            
            # perform the rebinning and create a new histo with the 
            # rebinned rates
            if histo.num_times > 0:
                edges = np.append(histo.tstart, histo.tstop[-1])
                new_counts = []
                for i in range(bin.num_chans):
                    new_cts, new_exposure, new_edges = method(
                        histo.counts[:, i],
                        histo.exposure,
                        edges, *args)
                    new_counts.append(new_cts)
                new_counts = np.array(new_counts).T

                new_histo = TimeEnergyBins(new_counts, new_edges[:-1],
                                           new_edges[1:], new_exposure, 
                                           bin.emin, bin.emax)
            else:
                new_histo = bin

            # now merge the split histo back together again
            histos_to_merge = [i for i in (pre, new_histo, post) if
                               i.num_times > 0]
            
            if len(histos_to_merge) > 0: 
                new_histos.append(TimeEnergyBins.merge_time(histos_to_merge))

        new_histo = TimeEnergyBins.merge_time(new_histos)

        return new_histo

    def slice_energy(self, emin, emax):
        """Perform a slice over an energy range and return a new TimeEnergyBins 
        object. Note that emin and emax values that fall inside a bin will 
        result in that bin being included.
        
        Args:
            emin (float): The start of the slice
            emax (float): The end of the slice
        
        Returns:           
            (:class:`TimeEnergyBins`)
        """
        mask = self._slice_energy_mask(emin, emax)
        obj = TimeEnergyBins(self.counts[:, mask], self.tstart, self.tstop,
                             self.exposure, self.emin[mask], self.emax[mask],
                             quality=self.quality)
        return obj

    def slice_time(self, tstart, tstop):
        """Perform a slice over a time range and return a new TimeEnergyBins 
        object. Note that tstart and tstop values that fall inside a bin will 
        result in that bin being included.
        
        Args:
            tstart (float): The start of the slice
            tstop (float): The end of the slice
        
        Returns:           
            (:class:`TimeEnergyBins`)
        """
        mask = self._slice_time_mask(tstart, tstop)
        cls = type(self)
        obj = cls(self.counts[mask, :], self.tstart[mask], self.tstop[mask], 
                  self.exposure[mask], self.emin, self.emax, 
                  quality=self.quality[mask])
        return obj

    @classmethod
    def merge_energy(cls, histos, **kwargs):
        """Merge multiple TimeEnergyBins together along the energy axis.
        
        Args:
            histos (list of :class:`TimeEnergyBins`): 
                A list containing the TimeEnergyBins to be merged
        
        Returns:
            (:class:`TimeEnergyBins`)
        """
        num = len(histos)
        # sort by channel edge
        emins = np.concatenate([[histo.emin[0]] for histo in histos])
        idx = np.argsort(emins)

        # concatenate the histos in order
        counts = histos[idx[0]].counts
        tstart = histos[idx[0]].tstart
        tstop = histos[idx[0]].tstop
        exposure = histos[idx[0]].exposure
        quality = histos[idx[0]].quality
        emin = histos[idx[0]].emin
        emax = histos[idx[0]].emax
        for i in range(1, num):
            bin_starts = histos[idx[i]].emin
            # make sure there is no overlap
            mask = (bin_starts >= emax[-1])
            if (~mask).sum() > 0:
                raise ValueError('Overlapping bins cannot be merged.  Only' \
                                 'non-overlapping bins can be merged.')

            counts = np.hstack((counts, histos[idx[i]].counts[:, mask]))
            emin = np.concatenate((emin, histos[idx[i]].emin[mask]))
            emax = np.concatenate((emax, histos[idx[i]].emax[mask]))

        # new TimeEnergyBins object
        merged_bins = cls(counts, tstart, tstop, exposure, emin, emax, 
                          quality=quality, **kwargs)
        return merged_bins

    @classmethod
    def merge_time(cls, histos, **kwargs):
        """Merge multiple TimeEnergyBins together along the time axis.
        
        Args:
            histos (list of :class:`TimeEnergyBins`): 
                A list containing the TimeEnergyBins to be merged
        
        Returns:
            (:class:`TimeEnergyBins`)
        """
        num = len(histos)
        # sort by start time
        tstarts = np.concatenate([[histo.tstart[0]] for histo in histos])
        idx = np.argsort(tstarts)

        # concatenate the histos in order
        counts = histos[idx[0]].counts
        tstart = histos[idx[0]].tstart
        tstop = histos[idx[0]].tstop
        exposure = histos[idx[0]].exposure
        emin = histos[idx[0]].emin
        emax = histos[idx[0]].emax
        quality = histos[idx[0]].quality
        for i in range(1, num):
            bin_starts = histos[idx[i]].tstart
            # make sure there is no overlap
            mask = (bin_starts >= tstop[-1])
            if (~mask).sum() > 0:
                raise ValueError('Overlapping bins cannot be merged.  Only' \
                                 'non-overlapping bins can be merged.')

            counts = np.vstack((counts, histos[idx[i]].counts[mask, :]))
            tstart = np.concatenate((tstart, histos[idx[i]].tstart[mask]))
            tstop = np.concatenate((tstop, histos[idx[i]].tstop[mask]))
            exposure = np.concatenate(
                (exposure, histos[idx[i]].exposure[mask]))
            quality = np.concatenate((quality, histos[idx[i]].quality[mask]))

        # new TimeEnergyBins object
        merged_bins = cls(counts, tstart, tstop, exposure, emin, emax, 
                          quality=quality, **kwargs)
        return merged_bins

    def _assert_range(self, valrange):
        assert valrange[0] <= valrange[1], \
            'Range must be in increasing order: (lo, hi)'
        return valrange
    
    def _calculate_good_segments(self, lo_edges, hi_edges):
        """Calculates the ranges of data that are contiguous segments
        
        Args:
            lo_edges (np.array): The lower bin edges
            hi_edges (np.array): The upper bin edges
        
        Returns:           
            ([(float, float), ...])
        """
        mask = (lo_edges[1:] != hi_edges[:-1])
        if mask.sum() == 0:
            return [(lo_edges[0], hi_edges[-1])]
        edges = np.concatenate(([lo_edges[0]], hi_edges[:-1][mask],
                                lo_edges[1:][mask], [hi_edges[-1]]))
        edges.sort()
        return edges.reshape(-1, 2).tolist()

    def _closest_edge(self, edges, val, which='either'):
        """Return the closest time bin edge
        
        Args:
            val (float): Input value
            which (str, optional): Options are: 
                
                * 'either' - closest edge to val; 
                * 'low' - closest edge lower than val; 
                * 'high' - closest edge higher than val. Default is 'either'
        
        Returns:           
            (float)
        """
        idx = np.argmin(np.abs(val - edges))
        if which == 'low':
            if (edges[idx] > val) and (idx - 1) >= 0:
                idx -= 1
        elif (which == 'high') and (idx + 1) < edges.size:
            if edges[idx] < val:
                idx += 1
        else:
            pass
        return edges[idx]

    def _slice_energy_mask(self, emin, emax):
        emin_snap = self.closest_energy_edge(emin, which='low')
        emax_snap = self.closest_energy_edge(emax, which='high')
        mask = (self.emin < emax_snap) & (self.emax > emin_snap)
        return mask

    def _slice_time_mask(self, tstart, tstop):
        tstart_snap = self.closest_time_edge(tstart, which='low')
        tstop_snap = self.closest_time_edge(tstop, which='high')
        mask = (self.tstart < tstop_snap) & (self.tstop > tstart_snap)
        return mask

    def __repr__(self):
        s = '<{0}: {1} time bins;\n'.format(self.__class__.__name__, 
                                            self.num_times)
        s += ' time range {0};\n'.format(self.time_range)
        if self._good_time_segments is not None:
            s += ' {0} time segments;\n'.format(len(self._good_time_segments))
        
        s += ' {0} energy bins;\n'.format(self.num_chans)
        s += ' energy range {0}'.format(self.energy_range)
        if self._good_energy_segments is not None:
            s += ';\n {0} energy segments'.format(len(self._good_energy_segments))
        
        return s+'>'


class ResponseMatrix():
    """A class defining a 2D detector response matrix, with an (input) photon
    axis and (output) channel axis.

    Parameters:
        matrix (np.array): The 2D matrix, of shape 
                           (num_photon_bins, num_channels)
        emin (np.array): The low edges of the (input) photon bins
        emax (np.array): The high edges of the (input) photon bins
        chanlo (np.array): The low edges of the (output) energy channels
        chanhi (np.array): The high edges of the (output) energy channels    
    """
    def __init__(self, matrix, emin, emax, chanlo, chanhi):
        try:
            iter(matrix)
            self._matrix = np.asarray(matrix)
        except:
            raise TypeError('matrix must be an iterable')
        if self._matrix.ndim != 2:
            raise TypeError('matrix must be a 2-dimensional array')
                
        try:
            iter(emin)
            self._emin = np.asarray(emin).flatten()
        except:
            raise TypeError('emin must be an iterable')
        try:
            iter(emax)
            self._emax = np.asarray(emax).flatten()
        except:
            raise TypeError('emax must be an iterable')
        
        try:
            iter(chanlo)
            self._chanlo = np.asarray(chanlo).flatten()
        except:
            raise TypeError('chanlo must be an iterable')
        try:
            iter(chanhi)
            self._chanhi = np.asarray(chanhi).flatten()
        except:
            raise TypeError('chanhi must be an iterable')

        if (self._emin.size != self._emax.size):
            raise ValueError('emin and emax must be the same length')        

        if (self._chanlo.size != self._chanhi.size):
            raise ValueError('chanlo and chanhi must be the same length')        
        
        if (self._matrix.shape[0] != self._emin.size) or \
           (self._matrix.shape[1] != self._chanlo.size):
            raise ValueError('matrix axis 0 must have same length as emin and '\
                             'emax.  matrix axis 1 must have same length as '\
                             'chanlo and chanhi.')
    
    @property
    def channel_centroids(self):
        """(np.array): The geometric mean of the energy channels"""
        return np.sqrt(self._chanlo * self._chanhi)

    @property
    def channel_widths(self):
        """(np.array): The energy channel widths"""
        return self._chanhi - self._chanlo

    @property
    def ebounds(self):
        """(:class:`Ebounds`): The energy bounds"""
        return Ebounds.from_bounds(self._chanlo, self._chanhi)

    @property
    def matrix(self):
        """(np.array): The raw matrix"""
        return self._matrix

    @property
    def num_chans(self):
        """(int): The number of energy channels"""
        return self._matrix.shape[1]

    @property
    def num_ebins(self):
        """(int): The number of photon bins"""
        return self._matrix.shape[0]
    
    @property
    def photon_bins(self):
        """(:class:`Ebounds`): The photon bins"""
        return Ebounds.from_bounds(self._emin, self._emax)

    @property
    def photon_bin_centroids(self):
        """(np.array): The geometric mean of the photon bins"""
        return np.sqrt(self._emin * self._emax)

    @property
    def photon_bin_widths(self):
        """(np.array): The photon bin widths"""
        return self._emax - self._emin

    def channel_effective_area(self):
        """Returns the effective area as a function of recorded channel energy
        (integrated over incident photon bins).
        
        Returns:
            (:class:`Bins`)
        """
        bins = Bins(self._matrix.sum(axis=0), self._chanlo, self._chanhi)
        return bins
    
    def effective_area(self, energy, interp_kind='linear'):
        """Calculate the effective area at a given energy or an array of
        energies.  This function interpolates the DRM (in log10(cm^2)) to
        provide the effective area for any energy within the photon energy
        bounds of the DRM.
        
        Args:
            energy (float or np.array): The photon energy or energies at which
                                        to calculate the effective area.
            interp_kind (str, optional): The kind of interpolation to be 
                                         passed to scipy.interp1d.  Default is
                                         'linear'.
        
        Returns:
            (np.array)
        """
        try:
            energy = np.asarray(energy, dtype=float)
        except:
            raise TypeError('energy must be a positive float')
                
        if np.any(energy < self.photon_bins.range[0]) or \
           np.any(energy > self.photon_bins.range[1]):
            raise ValueError('energy must be within the photon energy bounds ' \
                             'of the DRM') 
        
        effarea = self.photon_effective_area()
        diff_effarea = effarea.counts

        # create the DRM interpolator
        x = np.append(self._emin, self._emax[-1])
        y = np.append(diff_effarea, diff_effarea[-1])
        nz_mask = (y > 0)
        effarea_interpolator = interp1d(x[nz_mask], np.log10(y[nz_mask]), 
                                        kind=interp_kind, 
                                        fill_value='extrapolate')
        
        # do the interpolation
        effarea = 10.0**effarea_interpolator(energy)
        return effarea

    def fold_spectrum(self, function, params, channel_mask=None):
        """Fold a photon spectrum through a DRM to get a count spectrum
        
        Args: 
            function (<function>): 
                A photon spectrum function.  The function must accept a list of 
                function parameters as its first argument and an array of photon 
                energies as its second argument.  The function must return the 
                evaluation of the photon model in units of ph/s-cm^2-keV.
            params (list of float): A list of parameter values to be passed to
                                   the photon spectrum function
            channel_mask (np.array, optional): 
                A boolean mask where True indicates the channel is to be used 
                for folding and False indicates the channel is to not be used 
                for folding.  If omitted, all channels are used.
        
        Returns:        
            (np.array)
        """
        drm = self._matrix
        if channel_mask is not None:
            drm = drm[:, channel_mask]

        # evaluate photon model
        photon_model = function(params, self.photon_bin_centroids)

        # fold photon model through DRM
        counts = np.dot(drm.T, photon_model * self.photon_bin_widths)

        return counts
             
    def photon_effective_area(self):
        """Returns the effective area as a function of incident photon energy
        (integrated over recorded energy channels).
        
        Returns:        
            (:class:`Bins`)
        """
        bins = Bins(self._matrix.sum(axis=1), self._emin, self._emax)
        return bins

    def rebin(self, factor=None, edge_indices=None):
        """Rebins the channel energy axis of a DRM and returns a new response
        object.  Rebinning can only be used to downgrade the channel resolution
        and is constrained to the channel edges of the current DRM. 
        
        Rebinning can be performed by either defining an integral factor of the
        number of current energy channels to combine (e.g. factor=4 for 128
        energy channels would result in 32 energy channels), or by defining an
        index array into the current channel edges that will define the new set
        of edges.
        
        Args:
            factor (int, optional): The rebinning factor. Must set either this
                                    or `edge_indices`
            edge_indices (np.array, optional): The index array that represents
                                               which energy edges should remain
                                               in the rebinned DRM.
        
        Returns:
            (:class:`ResponseMatrix`)
        """
        if (factor is None) and (edge_indices is None):
            raise ValueError('Either factor or edge_indices must be set')
        
        elif factor is not None:
            try:
                factor = int(factor)
            except:
                raise TypeError('factor must be a positive integer')
            if factor < 1:
                raise ValueError('factor must be a positive integer')
            if (self.num_chans % factor) > 0:
                raise ValueError('factor must be factor of numchans: '\
                                 '{}'.format(self.num_chans))
            chanlo = self._chanlo.reshape(-1, factor)[:,0]
            chanhi = self._chanhi.reshape(-1, factor)[:,-1]
        
        elif edge_indices is not None:
            try:
                edge_indices = np.asarray(edge_indices)
            except:
                raise TypeError('new_edges must be a list, tuple or array')
            edge_indices = np.unique(edge_indices.astype(int))
            if (edge_indices[0] < 0) or (edge_indices[-1] > self.num_chans):
                raise ValueError('edge_indices outside valid range of ' \
                                 '0-{}'.format(self.num_chans))
            
            old_edges = np.append(self._chanlo, self._chanhi[-1])
            new_edges = old_edges[edge_indices]
            chanlo = new_edges[:-1]
            chanhi = new_edges[1:]
        
        # the new effective area matrix will just be summed over the channel
        # axis for the bins we are combining
        new_matrix = np.zeros((self.num_ebins, chanlo.size))
        sidx = (self._chanlo[:,np.newaxis] == chanlo[np.newaxis,:]).nonzero()[0]
        eidx = (self._chanhi[:,np.newaxis] == chanhi[np.newaxis,:]).nonzero()[0]
        for i in range(chanlo.size):
            new_matrix[:,i] = self._matrix[:,sidx[i]:eidx[i]+1].sum(axis=1)
        
        obj = type(self)(new_matrix, self._emin, self._emax, chanlo, chanhi)
        return obj   

    def resample(self, num_photon_bins=None, photon_bin_edges=None,
                 num_interp_points=20, interp_kind='linear'):
        """Resamples the incident photon axis of a DRM and returns a new 
        ResponseMatrix object.  Resampling can be used to downgrade the photon 
        energy resolution, upgrade the resolution, and/or redefine the edges of 
        the incident photon bins.  By definition, the resampling can only be 
        performed within the photon energy range of the current object.
        
        The process for resampling is to create a high-resolution grid for each
        new photon bin, interpolate the differential effective area onto that
        grid (interpolation is performed on log10(effective area)), and then
        integrate the differential effective area over the grid points in each
        bin.  The final effective area is then scaled by the ratio of the 
        initial number of photon bins to the requested number of photon bins.
        
        Args:
            num_photon_bins (int, optional): The number of photon bins in the
                                             new DRM. The bin edges will be 
                                             generated logarithmically. Only set
                                             this or `photon_bin_edges`.
            photon_bin_edges (np.array, optional): The array of photon bin edges.
                                                   Only set this or 
                                                   `num_photon_bins`
            num_interp_points (int, optional): The number of interpolation points
                                               used to integrate over for each
                                               new photon bin. Default is 20.
            interp_kind (str, optional): The kind of interpolation to be 
                                         passed to scipy.interp1d.  Default is
                                         'linear'.
        
        Returns:
            (:class:`ResponseMatrix`)
        """
        if (num_photon_bins is None) and (photon_bin_edges is None):
            raise ValueError('Either num_photon_bins or photon_bin_edges must' \
                             ' be set')
        elif num_photon_bins is not None:
            try:
                num_photon_bins = int(num_photon_bins)
                assert num_photon_bins > 0
            except:
                raise TypeError('num_photon_bins must be a positive integer')
            photon_bin_edges = np.geomspace(self._emin[0], self._emax[-1], 
                                            num_photon_bins+1)
        
        elif photon_bin_edges is not None:
            try:
                photon_bin_edges = np.asarray(photon_bin_edges)
            except:
                raise TypeError('photon_bin_edges must be an iterable')
            badmask = (photon_bin_edges < self._emin[0]) | \
                      (photon_bin_edges > self._emax[-1])
            if badmask.sum() > 0:
                raise ValueError('photon_bin_edges are beyond valid range of '\
                                 '{0:.2f}-{1:.2f}'.format(self._emin[0],
                                                          self._emax[-1]))
            photon_bin_edges = np.sort(photon_bin_edges)
            num_photon_bins = photon_bin_edges.size-1
        
        # differential effective area
        diff_effarea = self._matrix/self.photon_bin_widths[:,np.newaxis]
        
        # create an interpolation array over each new bin
        photon_bounds = np.vstack((photon_bin_edges[:-1], photon_bin_edges[1:]))
        interp_arr = np.geomspace(*photon_bounds, num_interp_points+2, axis=0)
        
        # create the DRM interpolator
        x = np.append(self._emin, self._emax[-1])
        y = np.vstack((diff_effarea, diff_effarea[-1,:]))
        y[y == 0] = 1e-10
        effarea_interpolator = interp1d(x, np.log10(y), axis=0, 
                                        kind=interp_kind, 
                                        fill_value='extrapolate')
        
        # interpolate the DRM over the photon interpolation array
        diff_effarea_interp = 10.**effarea_interpolator(interp_arr)
        badmask = np.isnan(diff_effarea_interp)
        diff_effarea_interp[badmask] = 0.0
        
        # integrate the DRM over the photon interpolation array
        diff_effarea_new = trapz(diff_effarea_interp, 
                                 interp_arr[:,:,np.newaxis], axis=0)
        
        # scale by ratio of photon bins
        drm = diff_effarea_new / (self.num_ebins/num_photon_bins)
        
        # create response object
        obj = type(self)(drm, photon_bin_edges[:-1], photon_bin_edges[1:], 
                         self._chanlo, self._chanhi)
        return obj   

    def __repr__(self):
        return '<ResponseMatrix: {0} energy bins; {1} ' \
               'channels>'.format(self.num_ebins, self.num_chans)


class Parameter():
    """A fit parameter class
    
    Parameters:
        value (float): The central fit value
        uncert (float or 2-tuple): The 1-sigma uncertainty. If a 2-tuple, then 
                                   is of the form (low, high)
        name (str, optional): The name of the parameter
        units (str, optional): The units of the parameter
        support (2-tuple, optional): The valid support of the parameter
    """
    def __init__(self, value, uncert, name='', units=None, 
                 support=(-np.inf, np.inf)):
        
        self._value = float(value)
        if isinstance(uncert, (tuple, list)):
            if len(uncert) == 2:
                pass
            elif len(uncert) == 1:
                uncert = (uncert[0], uncert[0])
            else:
                raise ValueError('uncertainty must be a 1- or 2-tuple')
        elif isinstance(uncert, float):
            uncert = (uncert, uncert)
        else:
            raise TypeError('uncertainty must be a float or 1- or 2-tuple')
        self._uncert = uncert
        
        self._units = units
        self._name = name
        self._support = support
    
    @property
    def name(self):
        """(str): The name of the parameter"""
        return self._name

    @property
    def support(self):
        """(2-tuple): The valid support of the parameter"""
        return self._support
    
    @property
    def uncertainty(self):
        """(float, float): The 1-sigma uncertainty"""
        return self._uncert
    
    @property
    def units(self):
        """(str): The units of the parameter"""
        return self._units

    @property
    def value(self):
        """(float): The central fit value"""
        return self._value

    def one_sigma_range(self):
        """Return the 1 sigma range of the parameter fit
        
        Returns:
            (tuple): 2-tuple (low, high)
        """
        return (self.value - self.uncertainty[0], 
                self.value + self.uncertainty[1])

    def to_fits_value(self):
        """Return as a tuple to be used for a FITS file

        Returns:
            (tuple): 3-value tuple (value, +uncertainty, -uncertainty)   
        """
        return (self.value, *self.uncertainty[::-1])

    def valid_value(self):
        """Check if the parameter value is within the allowed parameter range
        
        Returns:
            (bool)
        """
        if (self.value >= self.support[0]) and \
                (self.value <= self.support[1]):
            return True
        else:
            return False

    def _str_format(self):
        if (self.value > 0.005) and (self.uncertainty[0] > 0.005):
            value = '{0:.2f}'.format(self.value)
            uncertainty = tuple(
                ['{0:.2f}'.format(u) for u in self.uncertainty])
        else:
            value = '{0:.2e}'.format(self.value)
            val_coeff, val_exp = value.split('e')
            val_exp = int(val_exp)
            uncertainty = ['{0:.2e}'.format(u) for u in self.uncertainty]
            uncert_coeff = []
            uncert_exp = []
            for uncert in uncertainty:
                uncert_coeff.append(uncert.split('e')[0])
                uncert_exp.append(int(uncert.split('e')[1]))
        return (value, uncertainty)

    def __repr__(self):
        return '<{0}: {1}>'.format(self.__class__.__name__, self.name)
    
    def __str__(self):
        value, uncertainty = self._str_format()

        if uncertainty[0] == uncertainty[1]:
            s = '+/- {0}'.format(uncertainty[0])
        else:
            s = '+{1}/-{0}'.format(uncertainty[0], uncertainty[1])
        if self.units is None:
            return '{0}: {1} {2}'.format(self.name, value, s)
        else:
            return '{0}: {1} {2} {3}'.format(self.name, value, s, self.units)
