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

from .data_primitives import Ebounds, Gti, TimeEnergyBins, TimeChannelBins, \
                             EnergyBins
from .headers import FileHeaders
from .file import FitsFileContextManager
from .pha import Pha

__all__ = ['Phaii']

class Phaii(FitsFileContextManager):
    """PHAII class for time series of spectra.        
    """
    def __init__(self):
        super().__init__()
        self._ebounds = None
        self._data = None
        self._gti = None
        self._trigtime = None

    @property
    def data(self):
        """(:class:`~.data_primitives.TimeEnergyBins`): The PHAII data"""
        return self._data

    @property
    def ebounds(self):
        """(:class:`~.data_primitives.Ebounds`): The energy-channel mapping.
        If data does not have an ebounds, returns None."""
        if isinstance(self._data, TimeEnergyBins):
            return self._ebounds

    @property
    def energy_range(self):
        """(float, float): The energy range of the data.  If data does not have
        an ebounds, returns None."""
        if isinstance(self._data, TimeEnergyBins):
            return self._data.energy_range
    
    @property
    def gti(self):
        """(:class:`~.data_primitives.Gti`): The good time intervals"""
        return self._gti

    @property
    def num_chans(self):
        """(int): The number of energy channels"""
        return self._data.num_chans

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
        exposure = self._data.get_exposure(time_ranges=time_ranges)
        return exposure

    def rebin_energy(self, method, *args, energy_range=(None, None), **kwargs):
        """Rebin the PHAII in energy given a rebinning method. Produces a new 
        PHAII object.
        
        Note::
          If the data does not have an energy calibration (ebounds), then this
          function will bin by energy channels, and therefore the 
          ``energy_range`` argument should be a range of energy channels
          instead of energies.

        Args:
            method (<function>): The rebinning function
            *args: Arguments to be passed to the rebinning function
            energy_range ((float, float), optional): 
                The starting and ending energy (or channel) to rebin.  If 
                omitted, uses the full range of data.  Setting start or end to 
                ``None`` will use the data from the beginning or end of the data, respectively.
        Returns        
            (:class:`Phaii`)
        """
        emin, emax = self._assert_range(energy_range)
        if isinstance(self.data, TimeEnergyBins):
            data = self.data.rebin_energy(method, *args, emin=emin, emax=emax)
        else:
            data = self.data.rebin_channels(method, *args, chan_min=emin, 
                                            chan_max=emax)
        headers = self._build_headers(self.trigtime, *data.time_range, 
                                      data.num_chans)
        
        phaii = self.from_data(data, gti=self.gti, trigger_time=self.trigtime, 
                               headers=headers, **kwargs)
        return phaii

    def rebin_time(self, method, *args, time_range=(None, None), **kwargs):
        """Rebin the PHAII in time given a rebinning method. 
        Produces a new PHAII object.
        
        Args:
            method (<function>): The rebinning function
            *args: Arguments to be passed to the rebinning function
            time_range ((float, float), optional): 
                The starting and ending time to rebin.  If omitted, uses the 
                full range of data.  Setting start or end to ``None`` will use 
                the data from the beginning or end of the data, respectively.
        Returns        
            (:class:`Phaii`)
        """
        tstart, tstop = self._assert_range(time_range)
        data = self.data.rebin_time(method, *args, tstart=tstart, tstop=tstop)
        
        gti = Gti.from_list([data.time_range])
        gti = Gti.intersection(self.gti, gti)

        headers = self._build_headers(self.trigtime, *data.time_range, 
                                      data.num_chans)
        
        phaii = self.from_data(data, gti=self.gti, trigger_time=self.trigtime, 
                              headers=headers, **kwargs)
        return phaii
    
    def set_ebounds(self, ebounds):
        """Set the energy calibration (ebounds) of the data. If the data are
        not yet energy calibrated, this will convert the data object from
        :class:`~.data_primitives.TimeChannelBins` to 
        :class:`~.data_primitives.TimeEnergyBins``.  If the data already has an
        energy calibration, this method will update the calibration to with the
        new ebounds. The number of channels in ``ebounds`` must equal the number
        of channels of the data.
        
        Args:
            ebounds (:class:`~.data_primitives.Ebounds`): The ebounds
        """
        if not isinstance(ebounds, Ebounds):
            raise TypeError('ebounds must be an Ebounds object')
        
        if isinstance(self.data, TimeChannelBins):
            self._data = self.data.apply_ebounds(ebounds)
        else:
            if ebounds.num_intervals != self.num_chans:
                raise ValueError('Ebounds is of wrong size for this data')
            self._data = TimeEnergyBins(self.data.counts, self.data.tstart,
                                        self.data.tstop, self.data.exposure,
                                        ebounds.low_edges(), ebounds.high_edges())

        self._ebounds = ebounds
    
    def slice_energy(self, energy_ranges, **kwargs):
        """Slice the PHAII by one or more energy range. Produces a new PHAII 
        object.
        
        Note::
          If the data does not have an energy calibration (ebounds), then this
          function will slice by energy channels, and therefore the 
          ``energy_ranges`` argument should be a range(s) of energy channels
          instead of energies.

        Args:
            energy_ranges ([(float, float), ...]): 
                The energy ranges to slice the data to.
        
        Returns:        
            (:class:`Phaii`)
        """
        energy_ranges = self._assert_range_list(energy_ranges)
        
        if isinstance(self.data, TimeEnergyBins):
            data = [self.data.slice_energy(*self._assert_range(energy_range)) \
                    for energy_range in energy_ranges]
            data = TimeEnergyBins.merge_energy(data)
        else:
            data = [self.data.slice_channels(*self._assert_range(energy_range)) \
                    for energy_range in energy_ranges]
            data = TimeChannelBins.merge_channels(data)
        
        headers = self._build_headers(self.trigtime, *data.time_range, 
                                      data.num_chans)
        
        phaii = self.from_data(data, gti=self.gti, trigger_time=self.trigtime, 
                               headers=headers, **kwargs)
        
        return phaii

    def slice_time(self, time_ranges, **kwargs):
        """Slice the PHAII by one or more time range. Produces a new 
        PHAII object. The GTI will be automatically update to match the new
        time range(s).

        Args:
            time_ranges ([(float, float), ...]): 
                The time ranges to slice the data to.
        
        Returns:        
            (:class:`Phaii`)
        """
        time_ranges = self._assert_range_list(time_ranges)
        data = [self.data.slice_time(*self._assert_range(time_range)) \
                for time_range in time_ranges]
        
        gti = Gti.from_list(self.gti.as_list())
        for segment in data: 
            seg_gti = Gti.from_list([segment.time_range])
            gti = Gti.intersection(gti, seg_gti)
        
        if isinstance(self.data, TimeEnergyBins):       
            data = TimeEnergyBins.merge_time(data)
        else:
            data = TimeChannelBins.merge_time(data)
            
        headers = self._build_headers(self.trigtime, *data.time_range, 
                                      data.num_chans)
        
        phaii = self.from_data(data, gti=self.gti, trigger_time=self.trigtime, 
                               headers=headers, **kwargs)
        return phaii

    def to_lightcurve(self, time_range=None, energy_range=None,
                      channel_range=None):
        """Integrate the PHAII data over energy to produce a lightcurve.
        
        Note::
          If the data has not energy calibration, then ``energy_range`` is 
          ignored, and only ``channel_range`` is used.
        
        Args:
            time_range ((float, float), optional): 
                The time range of the lightcurve. If omitted, uses the entire 
                time range of the data.
            energy_range ((float, float), optional): 
                The energy range of the lightcurve. If omitted, uses the entire 
                energy range of the data.
            channel_range ((int, int), optional): 
                The channel range of the lightcurve. If omitted, uses the entire 
                energy range of the data.
        Returns:        
            (:class:`~.data_primitives.TimeBins`)
        """
        # slice to desired time range
        if time_range is not None:
            temp = self._data.slice_time(*self._assert_range(time_range))
        else:
            temp = self._data

        # limit integration to be over desired energy or channel range
        if channel_range is not None:
            self._assert_range(channel_range)
        
        if isinstance(self.data, TimeEnergyBins):
            if channel_range is not None:
                energy_range = (self._data.emin[channel_range[0]],
                                self._data.emax[channel_range[1]])
            if energy_range is not None:
                emin, emax = self._assert_range(energy_range)
            else:
                emin, emax = None, None

            # produce the TimeBins lightcurve object
            lc = temp.integrate_energy(emin=emin, emax=emax)
        
        else:
            if channel_range is None:
                channel_range = (None, None)
            lc = temp.integrate_channels(chan_min=channel_range[0],
                                         chan_max=channel_range[1])

        return lc
    
    def to_pha(self, time_ranges=None, energy_range=None, channel_range=None,
               **kwargs):
        """Integrate the PHAII data over time to produce a PHA object.
        
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
            **kwargs: Options passed to :meth:`.pha.Pha.from_data`
        
        Returns:        
            (:class:`~.pha.Pha`)
        """
        if isinstance(self.data, TimeChannelBins):
            raise RuntimeError('Energy calibration required to create a PHA object')
        
        if time_ranges is None:
            time_ranges = [self.time_range]
        time_ranges = self._assert_range_list(time_ranges)
        specs = []
        times = []
        for time_range in time_ranges:
            spec = self.to_spectrum(time_range=time_range,
                                    energy_range=energy_range,
                                    channel_range=channel_range)

            # PHA needs to have the full spectral range and zero out the ones
            # we're not using
            lomask = (self.data.emin < spec.range[0])
            pre_counts = np.zeros(np.sum(lomask), dtype=int)
            pre_lo_edges = self.data.emin[lomask]
            pre_hi_edges = self.data.emax[lomask]

            himask = (self.data.emax > spec.range[1])
            post_counts = np.zeros(np.sum(himask), dtype=int)
            post_lo_edges = self.data.emin[himask]
            post_hi_edges = self.data.emax[himask]

            counts = np.concatenate((pre_counts, spec.counts, post_counts))
            lo_edges = np.concatenate((pre_lo_edges, spec.lo_edges, 
                                       post_lo_edges))
            hi_edges = np.concatenate((pre_hi_edges, spec.hi_edges, 
                                       post_hi_edges))
            one_spec = EnergyBins(counts, lo_edges, hi_edges, spec.exposure[0])
            specs.append(one_spec)

            # need the actual time range
            times.append(self.data.slice_time(*time_range).time_range)
        
        # sum disjoint time selections
        data = EnergyBins.sum(specs)
        # values for the PHA object
        channel_mask = ~(lomask | himask)

        gti = Gti.from_list(times)

        pha = Pha.from_data(data, gti=gti, trigger_time=self.trigtime, 
                            channel_mask=channel_mask, **kwargs)
        
        return pha
    
    def to_spectrum(self, time_range=None, energy_range=None,
                    channel_range=None):
        """Integrate the PHAII data over time to produce a count spectrum

        Note::
          If the data has not energy calibration, then ``energy_range`` is 
          ignored, and only ``channel_range`` is used.

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
            
            if isinstance(self.data, TimeEnergyBins):
                if channel_range is not None:
                    energy_range = (self.data.emin[channel_range[0]],
                                    self.data.emax[channel_range[1]])
                temp = self.data.slice_energy(*self._assert_range(energy_range))
            else:
                if channel_range is None:
                    channel_range = self.data.channel_range
                temp = self.data.slice_channels(chan_min=channel_range[0],
                                                chan_max=channel_range[1])
        else:
            temp = self.data

        # limit integration to be over desired time range
        if time_range is not None:
            tstart, tstop = self._assert_range(time_range)
        else:
            tstart, tstop = None, None

        # produce the EnergyBins count spectrum
        spec = temp.integrate_time(tstart=tstart, tstop=tstop)
        return spec
    
    @classmethod
    def from_data(cls, data, gti=None, trigger_time=None, filename=None,
                  headers=None, **kwargs):
        """Create a PHAII object from a 
        :class:`~.data_primitives.TimeEnergyBins` data object.
        
        Args:
            data (:class:`~.data_primitives.TimeEnergyBins`): The PHAII data
            gti (:class:`~.data_primitives.Gti`, optional): 
                The Good Time Intervals object. If omitted, the GTI spans 
                (tstart, tstop) 
            trigger_time (float, optional): 
                The trigger time, if applicable. If provided, the data times 
                will be shifted relative to the trigger time.
            filename (str, optional): The name of the file
            headers (:class:`~.headers.FileHeaders`): The file headers
                 
        Returns:
            (:class:`Phaii`)
        """
        obj = cls()
        obj._filename = filename
        
        # set data and ebounds
        if not isinstance(data, TimeEnergyBins) and \
           not isinstance(data, TimeChannelBins):
            raise TypeError('data must be of type TimeEnergyBins')
        
        obj._data = data
        
        if isinstance(data, TimeEnergyBins):
            obj._ebounds = Ebounds.from_bounds(data.emin, data.emax)
        
        # set GTI
        if gti is not None:
            if not isinstance(gti, Gti):
                raise TypeError('gti must be of type Gti')
        else:
            gti = Gti.from_list([data.time_range])
        obj._gti = gti
                
        # update times to be relative to trigger time
        if trigger_time is not None:
            if trigger_time < 0.0:
                raise ValueError('trigger_time must be non-negative')
            obj._trigtime = trigger_time
        
        # set headers
        if headers is not None:
            if not isinstance(headers, FileHeaders):
                raise TypeError('headers must be of type FileHeaders')
        obj._headers = headers        
        
        return obj

    @classmethod
    def merge(cls, phaii_list):
        """Merge a list of Phaii objects into a new Phaii object.  The header
        from the first Phaii in the list is used in the new Phaii and 
        appropriately updated.
        
        Args:
            phaii_list (list): The list of Phaii objects to merge
                 
        Returns:
            (:class:`Phaii`)
        """
        # make sure all are Phaii objects and they all have the same number
        # of energy channels
        if any([not isinstance(phaii, Phaii) for phaii in phaii_list]):
            raise ValueError('phaii_list must be a list of Phaii objects')
        num_chans = phaii_list[0].num_chans
        for phaii in phaii_list[1:]:
            if phaii.num_chans != num_chans:
                raise ValueError('All Phaii objects must have the same number '\
                                 'of energy channels')
        
        # merge the data
        if isinstance(phaii_list[0].data, TimeEnergyBins):
            data = TimeEnergyBins.merge_time([phaii.data for phaii in phaii_list])
        else:
            data = TimeChannelBins.merge_time([phaii.data for phaii in phaii_list])
        
        # merge the GTIs
        gti = phaii_list[0].gti
        for phaii in phaii_list[1:]:
            gti = Gti.merge(gti, phaii.gti)
        
        # create header
        trigtime = phaii_list[0].trigtime
        headers = phaii_list[0]._build_headers(trigtime, *data.time_range, 
                                               data.num_chans)
        
        obj = cls.from_data(data, gti=gti, trigger_time=trigtime, 
                            filename=phaii_list[0].filename, headers=headers)
        return obj
    
    def _assert_exposure(self, exposure):
        """Sometimes the data files contain a negative exposure. That is very 
        very naughty, and we don't like dealing with naughty things.
        """
        mask = (exposure < 0.0)
        exposure[mask] = 0.0
        return exposure

    def _assert_range(self, valrange):
        if valrange[0] is None and valrange[1] is None:
            return valrange
        assert valrange[0] <= valrange[1], \
            'Range must be in increasing order: (lo, hi)'
        return valrange

    def _assert_range_list(self, range_list):
        try:
            iter(range_list[0])
        except:
            range_list = [tuple(range_list)]
        range_list = [self._assert_range(r) for r in range_list]
        return range_list

    def _build_headers(self, trigtime, tstart, tstop, num_chans):
        """This builds the headers for the FITS file.  This method needs
        to be specified in the inherited class.  The method should construct
        the headers from the minimum required arguments and additional keywords
        and return a :class:`~.headers.FileHeaders` object.
        
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
        
