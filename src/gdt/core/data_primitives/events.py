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
import copy
from .bins import ChannelBins, EnergyBins, TimeChannelBins, TimeEnergyBins
from .intervals import Ebounds

__all__ = ['EventList']


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

        self._events = np.empty(times.size, dtype=[('TIME', times.dtype.type),
                                                   ('PHA', channels.dtype.type)])
        self._events['TIME'] = times
        self._events['PHA'] = channels

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
        _, _, _, new_edges = method(counts, np.sqrt(counts), exposure, edges, *args)

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
