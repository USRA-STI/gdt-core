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

__all__ = ['BackgroundFitter']

import numpy as np
from gdt.core.data_primitives import EventList, TimeEnergyBins
from gdt.core.phaii import Phaii
from gdt.core.tte import PhotonList
from .primitives import BackgroundRates

class BackgroundFitter:
    """Class for fitting a background, given a fitting algorithm,
    to time-energy data (e.g. :class:`~gdt.core.phaii.Phaii`, 
    :class:`~gdt.core.tte.PhotonList`).
    
    When a BackgroundFitter is created, an algorithm must be specified.  In 
    particular, the algorithm must be a class that has two public methods:
    ``fit()`` and ``interpolate()``. 
    
    For PHAII data, the class, upon initialization, must take the following as 
    arguments: 
    
    * A 2D array of counts with shape (``num_times``, ``num_chans``);
    * A 1D array of time bin start times, shape (``num_times``,); 
    * A 1D array of time bin end times, shape (``num_times``,);
    * A 1D array of exposures for each time bin, shape (``num_times``,).

    While for TTE data, the class, upon initialization, must take the following 
    as an argument: 

    * A list, of length ``num_chans``, where each item is a numpy.ndarray 
      of event times.
    
    The ``fit()`` method takes no arguments, however, parameters that are 
    required for the algorithm may be specified as keywords.
    
    The interpolate() method must take in the following as arguments:
    
    * A 1D array of time bin start times, shape (``num_times``,);
    * A 1D array of time bin end times, shape (``num_times``,).
    
    Any additional parameters required can be specified as keywords. 
    The ``interpolate()`` method must return:
    
    * A 2D rates array, shape (``num_times``, ``num_chans``);
    * A 2D rate uncertainty array, shape (``num_times``, ``num_chans``).
    
    Additionally, the class can provide the following public attributes that 
    will be exposed by BackgroundFitter:
   
    * ``dof``: The degrees of freedom of the fits, array shape (``num_chans``,)
    * ``statistic``: The fit statistic for each fit, array shape (``num_chans``,)
    * ``statistic_name``: A string of the fit statistic used    
    """
    def __init__(self):
        self._data_obj = None
        self._method = None
        self._type = None
        self._statistic = None
        self._dof = None
        self._livetime = None
        self._parameters = None

    @property
    def dof(self):
        """(np.array): If available, the degrees-of-freedom of the fit for 
        each energy channel"""
        return getattr(self._method, 'dof', None)

    @property
    def livetime(self):
        """(float): The total livetime of the data used for the background"""
        return self._livetime

    @property
    def method(self):
        """(str): The name of the fitting algorithm class"""
        return self._method.__class__.__name__

    @property
    def parameters(self):
        """(dict): All parameters passed to the fitting algorithm"""
        return self._parameters

    @property
    def statistic(self):
        """(np.array): If available, the fit statistic for each energy 
        channel"""
        return getattr(self._method, 'statistic', None)

    @property
    def statistic_name(self):
        """(str): If available, the name of the fit statistic"""
        return getattr(self._method, 'statistic_name', None)

    @property
    def type(self):
        """(str): The type of background algorithm, either 'binned' or 
        'unbinned'"""
        return self._type

    def fit(self, **kwargs):
        """Perform a background fit of the data
        
        Args:
            **kwargs: Options to be passed as parameters to the fitting class
        """
        self._parameters = kwargs
        self._method.fit(**kwargs)

    def interpolate_bins(self, tstart, tstop, **kwargs):
        """Interpolate the fitted background model over a set of bins.
        The exposure is calculated for each bin of the background model 
        in case the background model counts is needed.
        
        Args:
            tstart (np.array): The starting times
            tstop (np.array): The ending times
            **kwargs: Options to be passed as parameters to the interpolation 
                      method 
        
        Returns:
            (:class:`~gdt.background.primitives.BackgroundRates`)
        """
        # do the interpolation
        rate, rate_uncert = self._method.interpolate(tstart, tstop, **kwargs)
        # get the exposure
        numtimes = tstart.shape[0]
        exposure = np.array([self._data_obj.get_exposure((tstart[i], tstop[i])) \
                             for i in range(numtimes)])
        # create the rates object
        rates = BackgroundRates(rate, rate_uncert, tstart, tstop,
                                self._data_obj.data.emin,
                                self._data_obj.data.emax, exposure=exposure)

        return rates

    def interpolate_times(self, times, **kwargs):
        """Interpolate the fitted background model over a set of times.
        Does not calculate an exposure since this returns a set of point
        estimates of the background rates.
        
        Args:
            tstart (np.array): The sampling times
            **kwargs: Options to be passed as parameters to the interpolation 
                      method 
        
        Returns:
            (:class:`~gdt.background.primitives.BackgroundRates`)
        """
        # do the interpolation
        rate, rate_uncert = self._method.interpolate(times, times, **kwargs)

        # create the rates object
        rates = BackgroundRates(rate, rate_uncert, times, times,
                                self._data_obj.ebounds.low_edges(),
                                self._data_obj.ebounds.high_edges())

        return rates

    @classmethod
    def from_phaii(cls, phaii, method, time_ranges=None):
        """Create a background fitter from a PHAII object
        
        Args:
            phaii (:class:`~gdt.core.phaii.Phaii`): A PHAII data object
            method (<class>): A background fitting/estimation class for binned data 
            time_ranges ([(float, float), ...]): 
                The time range or time ranges over which to fit the background. 
                If omitted, uses the full time range of the data 
        
        Returns:        
            (:class:`BackgroundFitter`)
        """
        if not isinstance(phaii, Phaii):
            raise TypeError('Input data must be a Phaii object')

        obj = cls()
        obj._data_obj = phaii
        obj._validate_method(method)
        time_ranges = obj._validate_time_ranges(time_ranges)

        # Slice the PHAII data and merge if multiple slices
        data = [phaii.data.slice_time(trange[0], trange[1]) for trange in
                time_ranges]
        data = TimeEnergyBins.merge_time(data)
        obj._method = method(data.counts, data.tstart, data.tstop,
                             data.exposure)
        obj._type = 'binned'
        obj._livetime = np.sum(data.exposure)
        return obj

    @classmethod
    def from_tte(cls, tte, method, time_ranges=None):
        """Create a background fitter from a TTE object
        
        Args:
            tte (:class:`~gdt.core.tte.PhotonList`): A PhotonList data object
            method (<class>): A background fitting/estimation class for unbinned data 
            time_ranges ([(float, float), ...]): 
                The time range or time ranges over which to fit the background. 
                If omitted, uses the full time range of the data 
        
        Returns:        
            (:class:`BackgroundFitter`)
        """
        if not isinstance(tte, PhotonList):
            raise TypeError('Input data must be a PhotonList object')

        obj = cls()
        obj._data_obj = tte
        obj._validate_method(method)
        time_ranges = obj._validate_time_ranges(time_ranges)

        # Slice the TTE data and merge if multiple slices
        data = [tte.data.time_slice(trange[0], trange[1]) for trange in
                time_ranges]
        data = EventList.merge(data)
        data.sort_time()
        # pull out the events in each channel
        events = [data.channel_slice(i, i).times for i in
                  range(tte.num_chans)]

        obj._method = method(events)
        obj._type = 'unbinned'
        obj._livetime = data.get_exposure(time_ranges=time_ranges)
        return obj

    def _validate_method(self, method):
        try:
            method
        except:
            raise NameError('Input method is not a known function')
        
        try:
            has_fit = callable(method.fit)
            has_interp = callable(method.interpolate)
        except:
            raise NotImplementedError(
                "User-defined Background class must have "
                "both fit() and an interpolate() methods")

        if (not has_fit) or (not has_interp):
            raise NotImplementedError(
                "User-defined Background class must have "
                "both fit() and an interpolate() methods")

    def _validate_time_ranges(self, time_ranges):
        if time_ranges is None:
            time_ranges = [self._data_obj.time_range]
        try:
            iter(time_ranges[0])
        except:
            raise TypeError('time_ranges must be a list of tuples')
        return time_ranges

    def __repr__(self):
        s = "<BackgroundFitter: {0};\n".format(self.method)
        s += '{} data;'.format(self._type)
        try:
            s += '\n{0}/dof: {1:.2f}/{2}'.format(self.statistic_name, 
                                                 self.statistic.sum(),
                                                 int(self.dof.sum()))
        except:
            pass                              
        return s + '>'
        
