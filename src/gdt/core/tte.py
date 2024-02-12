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
import os
import numpy as np

from .data_primitives import Ebounds, Gti, EventList
from .file import FitsFileContextManager
from .headers import FileHeaders
from .phaii import Phaii
from .pha import Pha

__all__ = ['PhotonList']

class PhotonList(FitsFileContextManager):
    """Class for Photon Lists or Time-Tagged Event data
    """
    def __init__(self):
        super().__init__()
        self._data = None
        self._gti = None
        self._trigtime = None
        self._event_deadtime = 0.0
        self._overflow_deadtime = 0.0
    
    @property
    def data(self):
        """(:class:`~.data_primitives.EventList`): The event data"""
        return self._data
    
    @property
    def ebounds(self):
        """(:class:`~.data_primitives.Ebounds`): The energy-channel mapping"""
        return self.data.ebounds

    @property
    def energy_range(self):
        """(float, float): The energy range of the data"""
        return self._data.energy_range
    
    @property
    def event_deadtime(self):
        """(float): The deadtime imposed per event"""
        return self._event_deadtime
    
    @property
    def gti(self):
        """(:class:`~.data_primitives.Gti`): The good time intervals"""
        return self._gti
    
    @property
    def num_chans(self):
        """(int): The number of energy channels"""
        return self._data.num_chans

    @property
    def overflow_deadtime(self):
        """(float): The per-event deadtime imposed by events in the overflow 
        channel"""
        return self._overflow_deadtime
    
    @property
    def time_range(self):
        """(float, float): The time range of the data"""
        return self._data.time_range

    @property
    def trigtime(self):
        """(float): The trigger time of the data, if available."""
        return self._trigtime

    def get_exposure(self, time_ranges=None):
        """Calculate the total exposure of a time range or time ranges of data.

        Args:
            time_ranges ([(float, float), ...], optional): 
                The time range or time ranges over which to calculate the 
                exposure. If omitted, calculates the total exposure of the data.
        
        Returns:        
            (float)
        """
        if time_ranges is None:
            time_ranges = self.time_range
        time_ranges = self._assert_range_list(time_ranges)
        exposure = self._data.get_exposure(time_ranges=time_ranges, 
                                           event_deadtime=self.event_deadtime,
                                           overflow_deadtime=self.overflow_deadtime)
        return exposure

    def rebin_energy(self, method, *args, **kwargs):
        """Rebin the PhotonList in energy given a rebinning method. 
        Produces a new PhotonList object.

        Args:
            method (<function>): The rebinning function
            *args: Arguments to be passed to the rebinning function
        
        Returns        
            (:class:`PhotonList`)
        """
        data = self._data.rebin_energy(method, *args, 
                                       event_deadtime=self.event_deadtime,
                                       overflow_deadtime=self.overflow_deadtime)
        headers = self._build_headers(self.trigtime, *data.time_range, 
                                      data.num_chans)
        obj = self.from_data(data, gti=self.gti, trigger_time=self.trigtime, 
                             headers=headers, 
                             event_deadtime=self.event_deadtime,
                             overflow_deadtime=self.overflow_deadtime, **kwargs)
        return obj

    def set_ebounds(self, ebounds):
        """Set the energy calibration (ebounds) of the data. If the data already 
        has an energy calibration, this method will update the calibration to 
        the new ebounds.
        
        Args:
            ebounds (:class:`~.data_primitives.Ebounds`): The ebounds
        """
        self.data.ebounds = ebounds        

    def slice_energy(self, energy_ranges, **kwargs):
        """Slice the PhotonList by one or more energy ranges. Produces a new 
        PhotonList object.

        Note::
          If the data does not have an energy calibration (ebounds), then this
          function will slice by energy channels, and therefore the 
          ``energy_ranges`` argument should be a range(s) of energy channels
          instead of energies.
        
        Args:
            energy_ranges ([(float, float), ...]): 
                The energy ranges to slice the data to.

        Returns:
            (:class:`PhotonList`)
        """
        energy_ranges = self._assert_range_list(energy_ranges)
        if self.ebounds is not None:
            data = [self.data.energy_slice(*self._assert_range(energy_range)) \
                    for energy_range in energy_ranges]
        else:
            data = [self.data.channel_slice(*self._assert_range(energy_range)) \
                    for energy_range in energy_ranges]
        data = EventList.merge(data, sort=True)

        headers = self._build_headers(self.trigtime, *data.time_range, 
                                      data.num_chans)
        tte = self.from_data(data, gti=self.gti, trigger_time=self.trigtime,
                             headers=headers, 
                             event_deadtime=self.event_deadtime,
                             overflow_deadtime=self.overflow_deadtime, **kwargs)
        return tte

    def slice_time(self, time_ranges, **kwargs):
        """Slice the PhotonList by one or more time range. Produces a new 
        PhotonList object.
        
        Args:
            time_ranges ([(float, float), ...]): 
                The time ranges to slice the data to.
        
        Returns:
            (:class:`PhotonList`)
        """
        time_ranges = self._assert_range_list(time_ranges)
        data = [self.data.time_slice(*self._assert_range(time_range)) \
                for time_range in time_ranges]
        
        gti = Gti.from_list(self.gti.as_list())
        for segment in data: 
            seg_gti = Gti.from_list([segment.time_range])
            gti = Gti.intersection(gti, seg_gti)
        
        data = EventList.merge(data, sort=True)

        headers = self._build_headers(self.trigtime, *data.time_range, 
                                      data.num_chans)        
        tte = self.from_data(data, gti=gti, trigger_time=self.trigtime,
                             headers=headers, 
                             event_deadtime=self.event_deadtime,
                             overflow_deadtime=self.overflow_deadtime, **kwargs)
        
        return tte
    
    def to_pha(self, time_ranges=None, energy_range=None, channel_range=None,
               **kwargs):
        """Integrate the PhotonList data over one or more time ranges to 
        produce a PHA object

        Note::
          If the data does not have an energy calibration (ebounds), then a 
          PHA object cannot be created and calling this method will raise an
          exception.

        Args:
            time_ranges ([(float, float), ...], optional): 
                The time range of the spectrum. If omitted, uses the entire 
                time range of the data.
            energy_range ((float, float), optional): 
                The energy range of the spectrum. If omitted, uses the entire 
                energy range of the data.
            channel_range ((int, int), optional): 
                The channel range of the spectrum. If omitted, uses the entire 
                energy range of the data.

            **kwargs: Options passed to :meth:`.pha.PHA.from_data`
        
        Returns:
            (:class:`~.pha.PHA`)
        """
        if self.ebounds is None:
            raise RuntimeError('Energy calibration required to create a PHA object')

        if time_ranges is None:
            time_ranges = [self.time_range]
        time_ranges = self._assert_range_list(time_ranges)

        segs = [self._data.time_slice(*self._assert_range(time_range)) \
                for time_range in time_ranges]
        times = [seg.time_range for seg in segs]
        el = EventList.merge(segs)
        if energy_range is not None:
            el = el.energy_slice(*self._assert_range(energy_range))
            lomask = (np.array(self.ebounds.high_edges()) < energy_range[0])
            himask = (np.array(self.ebounds.low_edges()) > energy_range[1])
        elif channel_range is not None:
            el = el.channel_slice(*self._assert_range(channel_range))
            idx = np.arange(self.num_chans)
            lomask = (idx < channel_range[0])
            himask = (idx > channel_range[1])
        else:
            lomask = himask = np.zeros(self.num_chans, dtype=bool)
        
        data = el.count_spectrum(event_deadtime=self.event_deadtime,
                                 overflow_deadtime=self.overflow_deadtime)
        # values for the PHA object
        channel_mask = ~(lomask | himask)
        
        gti = Gti.from_list(times)
        
        pha = Pha.from_data(data, gti=gti, trigger_time=self.trigtime, 
                            channel_mask=channel_mask, **kwargs)

        return pha

    def to_phaii(self, bin_method, *args, time_range=None, energy_range=None,
                 channel_range=None, phaii_class=Phaii, headers=None, **kwargs):
        """Convert the PhotonList data to PHAII data by binning the data in 
        time.

        Note::
          If the data has no energy calibration, then ``energy_range`` is 
          ignored, and only ``channel_range`` is used.
        
        Args:
            bin_method (<function>): A binning function for unbinned data
            *args: Arguments to pass to the binning function
            time_range ([(float, float), ...], optional):
                The time range of the spectrum. If omitted, uses the entire 
                time range of the data.
            energy_range ((float, float), optional): 
                The energy range of the spectrum. If omitted, uses the entire 
                energy range of the data.
            channel_range ((int, int), optional): 
                The channel range of the spectrum. If omitted, uses the entire 
                energy range of the data.
            phaii_class (class): The Phaii subclass that the data will be 
                                 converted to.  Default is the base 
                                 :class:`~.phaii.Phaii` class.
            headers (:class:`~.headers.FileHeaders`, optional): 
                The PHAII headers 
            **kwargs: Options to pass to the binning function
        
        Returns:
            (:class:`~.phaii.Phaii`)
        """
        if not issubclass(phaii_class, Phaii):
            raise TypeError('phaii_class must be set to a functional Phaii ' \
                            'derived class.')
        
        # slice to desired energy or channel range
        if (channel_range is not None) or (energy_range is not None):
            if channel_range is not None:
                self._assert_range(channel_range)

            if self.ebounds is not None:
                if channel_range is not None:
                    energy_range = (self.ebounds.low_edges()[channel_range[0]],
                                    self.ebounds.high_edges()[channel_range[1]])
            else:
                if channel_range is None:
                    channel_range = self.data.channel_range
                energy_range = channel_range
                    
            obj = self.slice_energy(energy_ranges=self._assert_range(energy_range))
        else:
            obj = self
        
        # slice to desired time range
        if time_range is None:
            pass
        else:
            obj = obj.slice_time(time_range)

        # do the time binning to create the TimeEnergyBins or TimeChannelBins
        bins = obj.data.bin(bin_method, *args, **kwargs)
        if (energy_range is not None) and (self.ebounds is not None):
            bins = bins.slice_energy(*energy_range)
        
        phaii = phaii_class.from_data(bins, gti=obj.gti, 
                                      trigger_time=obj.trigtime,
                                      headers=headers, **kwargs)
        return phaii

    def to_spectrum(self, time_range=None, energy_range=None,
                    channel_range=None):
        """Integrate the PhotonList data over time to produce a count spectrum

        Args:
            time_range ((float, float), optional):
                The time range of the spectrum. If omitted, uses the entire 
                time range of the data.
            energy_range ((float, float), optional): 
                The energy range of the spectrum. If omitted, uses the entire 
                energy range of the data.
            channel_range ((int, int), optional): 
                The channel range of the spectrum. If omitted, uses the entire 
                energy range of the data.

        Returns:
            (:class:`~.data_primitives.EnergyBins`)
        """
        # slice to desired energy or channel range   
        if (channel_range is not None) or (energy_range is not None):
            if channel_range is not None:
                self._assert_range(channel_range)
                
            if self.ebounds is not None:
                if channel_range is not None:
                    energy_range = (self.ebounds.low_edges()[channel_range[0]],
                                    self.ebounds.high_edges()[channel_range[1]])
                temp = self.data.energy_slice(*self._assert_range(energy_range))
            else:
                if channel_range is None:
                    channel_range = self.data.channel_range
                temp = self.data.channel_slice(*channel_range)
        else:
            temp = self._data
        
        # time slice
        if time_range is None:
            time_range = self.time_range
        segment = temp.time_slice(*self._assert_range(time_range))
        
        spec = segment.count_spectrum()
        return spec
       
    @classmethod
    def from_data(cls, data, gti=None, trigger_time=None, filename=None,
                  headers=None, event_deadtime=0.0, overflow_deadtime=0.0,
                  **kwargs):
        """Create a PhotonList object from an 
        :class:`~.data_primitives.EventList` object.

        Args:
            data (:class:`~.data_primitives.EventList`): The event data
            gti (:class:`~.data_primitives.Gti`, optional): 
                The Good Time Intervals object. If omitted, the GTI spans 
                (tstart, tstop) 
            trigger_time (float, optional): 
                The trigger time, if applicable. If provided, the data times 
                will be shifted relative to the trigger time.
            filename (str): The filename
            headers (:class:`~.headers.FileHeaders`): The file headers
            event_deadtime (float, optional): The deadtime imposed per event.
                                              Default is 0.0.
            overflow_deadtime (float, optional): The per-event deadtime imposed 
                                                 by events in the overflow 
                                                 channel. Default is 0.0.
        
        Returns:
            (:class:`PhotonList`)
        """
        obj = cls()
        obj._filename = filename
        
        # set data
        if not isinstance(data, EventList):
            raise TypeError('data must be of type EventList')
        obj._data = data
        
        # set GTI
        if gti is not None:
            if not isinstance(gti, Gti):
                raise TypeError('gti must be of type Gti')
        else:
            gti = Gti.from_list([data.time_range])
        obj._gti = gti
                
        # set headers
        if headers is not None:
            if not isinstance(headers, FileHeaders):
                raise TypeError('headers must be of type FileHeaders')
        obj._headers = headers        
 
        if trigger_time is not None:
            if trigger_time < 0.0:
                raise ValueError('trigger_time must be non-negative')
            obj._trigtime = trigger_time
        
        try:
            event_deadtime = float(event_deadtime)
            overflow_deadtime = float(overflow_deadtime)
        except:
            raise TypeError('deadtime must be a positive float')
        if event_deadtime < 0.0 or overflow_deadtime < 0.0:
            raise ValueError('deadtime must be non-negative')
        obj._event_deadtime = event_deadtime
        obj._overflow_deadtime = overflow_deadtime
       
        return obj

    @classmethod
    def merge(cls, ttes, primary=0, force_unique=True):
        """Merge a list of PhotonList objects into a single new PhotonList o
        bject.
        
        Warning: 
            The amount of data in a single PhotonList object can be quite large. 
            It is up to you to determine if you have enough memory to support 
            the merge.
        
        Args:
            ttes (list): A list of PhotonList objects to merge
            primary (int): 
                The index into the list of PhotonList objects to designate the 
                primary PhotonList object.  The primary object will be the 
                reference for the header information for the new merged 
                PhotonList object. Default is the first PhotonList object in 
                the list (primary=0).
            force_unique (bool, optional): 
                If True, force all events to be unique via brute force sorting. 
                If False, the EventLists will only be checked and masked for 
                overlapping time ranges. Events can potentially be lost if the 
                merged EventLists contain overlapping times (but not necessarily 
                duplicate events), however this method is much faster.  
                Default is True.
        
        Returns:
            (:class:`PhotonList`)
        """
        # make sure all are PhotonList objects and they all have the same number
        # of energy channels
        if any([not isinstance(tte, PhotonList) for tte in ttes]):
            raise ValueError('ttes must be a list of PhotonList objects')
        num_chans = ttes[0].num_chans
        for tte in ttes[1:]:
            if tte.num_chans != num_chans:
                raise ValueError('All PhotonList objects must have the same ' \
                                  'number of energy channels')

        if primary < 0 or primary > len(ttes)-1:
            raise ValueError('primary out of range for {} TTE objects' \
                             ''.format(len(ttes)))
        
        # merge the Event data
        data = [tte.data for tte in ttes]
        data = EventList.merge(data, sort=True, force_unique=force_unique)

        # merge the GTIs
        gtis = [tte.gti for tte in ttes]
        merged_gti = gtis[0]
        for gti in gtis[1:]:
            merged_gti = Gti.merge(merged_gti, gti)
            
        trigtime = ttes[primary].trigtime

        obj = cls.from_data(data, gti=merged_gti, trigger_time=trigtime,
                            headers=ttes[primary].headers, 
                            event_deadtime=ttes[primary].event_deadtime,
                            overflow_deadtime=ttes[primary].overflow_deadtime)
        return obj

    def _assert_range(self, valrange):
        assert valrange[0] <= valrange[1], \
            'Range must be in increasing order: (lo, hi)'
        return valrange

    def _assert_range_list(self, range_list):
        try:
            iter(range_list[0])
        except:
            range_list = [range_list]
        range_list = [self._assert_range(r) for r in range_list]
        return range_list

    def _build_headers(self, trigtime, tstart, tstop, num_chans):
        """This builds the headers for the FITS file.  This method needs
        to be specified in the inherited class.  The method should construct
        the headers from the minimum required arguments and additional keywords
        and return a :class:`FileHeaders` object.
        
        Args:
            trigtime (float or None): The trigger time.  Set to None if no
                                          trigger time.
            tstart (float): The start time
            tstop (float): The stop time
            num_chans (int): Number of detector energy channels
        
        Returns:
            (:class:`~.headers.FileHeaders`)
        """
        pass

    def __repr__(self):
        s = '<{0}: '.format(self.__class__.__name__)
        if self.filename is not None:
            s += '{};'.format(self.filename)
        if self.trigtime is not None:
            s += '\n trigger time: {};'.format(self.trigtime)
        s += '\n time range {};'.format(self.time_range)
        s += '\n energy range {}>'.format(self.energy_range)
        return s
