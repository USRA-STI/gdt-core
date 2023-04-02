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
from gdt.core.background.primitives import BackgroundSpectrum
from gdt.core.tte import PhotonList
from gdt.core.data_primitives import EventList, Ebounds, Gti
from gdt.core.response import Rsp
from .generators import *

__all__ = ['TteBackgroundSimulator', 'TteSourceSimulator']


class TteSourceSimulator:
    """Simulate TTE or EventList data for a source spectrum given a detector 
    response, spectral model and time profile model.  The spectral shape is 
    fixed throughout the time profile of the signal, but the amplitude of the 
    spectrum is time-dependent, set by the time profile function.

    Parameters:
        rsp (:class:`~gdt.core.response.Rsp`): A detector response object
        spec_func (:class:`~gdt.spectra.functions.Function`):
            A photon model function
        spec_params (iterable): The parameters for the function
        time_func (<function>): A time profile function
        time_params (iterable): Parameters for the time profile function
        sample_period (float, optional): The sampling period of the simulator
            in seconds. Default is 0.001. The simulator will produce arrival 
            times consistent with a spectrum over a finite time slice.  This 
            time slice should be short enough to approximate a non-homogeneous 
            Poisson process, but long enough to allow for a tractable 
            computation time.
        deadtime (float, optional): The dead time in seconds for each recorded 
                                    count during which another count cannot be 
                                    recorded. Default 0.
    """
    def __init__(self, rsp, spec_func, spec_params, time_func, time_params,
                 sample_period=0.001, deadtime=0.0):
        
        
        if sample_period <= 0.0:
            raise ValueError('Sample period must be positive')
        if deadtime < 0.0:
            raise ValueError('Deadtime must be non-negative')
        
        self._rsp = rsp
        self._spec_func = spec_func
        self._spec_params = spec_params
        self._time_func = time_func
        self._time_params = time_params
        self._sample_per = sample_period
        self._spec_gen = VariableSourceSpectrumGenerator(rsp, spec_func,
                                                         spec_params,
                                                         sample_period)
        self._event_gen = EventSpectrumGenerator(np.zeros(rsp.num_chans),
                                                 self._sample_per, 
                                                 min_sep=deadtime)

    def set_response(self, rsp):
        """Set/change the detector response.
        
        Args:
            rsp (:class:`~gdt.core.response.RSP`): A detector response object        
        """
        if not isinstance(rsp, Rsp):
            raise TypeError('rsp must be a Rsp object')
        self._rsp = rsp
        self._spec_gen = VariableSourceSpectrumGenerator(rsp, self._spec_func,
                                                         self._spec_params,
                                                         self._sample_per)

    def set_spectrum(self, spec_func, spec_params):
        """Set/change the spectrum.
        
        Args:
            spec_func (:class:`~gdt.spectra.functions.Function`):
                A photon model function
            spec_params (iterable): The parameters for the function
        """
        self._spec_func = spec_func
        self._spec_params = spec_params
        self._spec_gen = VariableSourceSpectrumGenerator(self._rsp,
                                                         self._spec_func,
                                                         self._spec_params,
                                                         self._sample_per)

    def set_time_profile(self, time_func, time_params):
        """Set/change the time profile.
        
        Args:
            time_func (<function>): A time profile function
            time_params (iterable): Parameters for the time profile function
        """
        self._time_func = time_func
        self._time_params = time_params

    def simulate(self, tstart, tstop):
        """Generate an EventList containing the individual counts from the 
        simulation.
        
        Args:
            tstart (float): The start time of the simulation
            tstop (float): The stop time of the simulation
        
        Returns:
            (:class:`~gdt.core.data_primitives.EventList`)
        """
        # create the time grid
        dur = (tstop - tstart)
        numpts = int(round(dur / self._sample_per))
        time_array = np.linspace(tstart, tstop, numpts)

        # calculate the spectral amplitudes over the grid
        amps = self._time_func(time_array, *self._time_params)

        times = []
        chans = []
        for i in range(numpts):
            # update amplitude and generate the count spectrum
            self._spec_gen.amp = amps[i]
            self._event_gen.spectrum = next(self._spec_gen).counts
            # generate the count arrival times for the time slice spectrum
            events = next(self._event_gen)
            if events is not None:
                times.extend((events[0] + time_array[i]).tolist())
                chans.extend(events[1].tolist())

        # create the eventlist
        eventlist = EventList(times=times, channels=chans, 
                              ebounds=self._rsp.ebounds)
        eventlist.sort_time()
        return eventlist

    def to_tte(self, tstart, tstop, **kwargs):
        """Generate a TTE object containing the individual counts from the 
        simulation.
        
        Args:
            tstart (float): The start time of the simulation
            tstop (float): The stop time of the simulation
            **kwargs: Options to pass to :class:`~gdt.core.tte.PhotonList`
        
        Returns:
            (:class:`~gdt.core.tte.PhotonList`)
        """
        eventlist = self.simulate(tstart, tstop)
        gti = Gti.from_bounds([eventlist.time_range[0]], 
                              [eventlist.time_range[1]])
        tte = PhotonList.from_data(eventlist, gti=gti, **kwargs)
        return tte


class TteBackgroundSimulator:
    """Simulate TTE or EventList data given a modeled background spectrum and
    time profile model. The spectrum is fixed throughout the time profile of
    the background, but the amplitude of the background is time-dependent, set
    by the time profile function.

    Parameters:
        bkgd_spectrum (:class:`~gdt.background.primitives.BackgroundSpectrum`):
            A modeled background spectrum
        distrib (str): The distribution from which the background is
                       simulated; either 'Poisson' or 'Gaussian'
        time_func (<function>): A time profile function
        time_params (iterable): Parameters for the time profile function
        sample_period (float, optional): The sampling period of the simulator
            in seconds. Default is 0.001. The simulator will produce arrival 
            times consistent with a spectrum over a finite time slice.  This 
            time slice should be short enough to approximate a non-homogeneous 
            Poisson process, but long enough to allow for a tractable 
            computation time.
        deadtime (float, optional): The dead time in seconds for each recorded 
                                    count during which another count cannot be 
                                    recorded. Default is 0.
    """
    def __init__(self, bkgd_spectrum, distrib, time_func, time_params,
                 sample_period=0.001, deadtime=0.0):
        
        if sample_period <= 0.0:
            raise ValueError('Sample period must be positive')
        if deadtime < 0.0:
            raise ValueError('Deadtime must be non-negative')
        
        self._spec_gen = None
        self._bkgd = bkgd_spectrum
        self._time_func = time_func
        self._time_params = time_params
        self._sample_per = sample_period
        self._deadtime = deadtime
        self._event_gen = EventSpectrumGenerator(np.zeros(self._bkgd.size),
                                                 self._sample_per, 
                                                 min_sep=deadtime)
        self.set_background(bkgd_spectrum, distrib)

    def set_background(self, bkgd_spectrum, distrib):
        """Set/change the spectrum.
        
        Args:
            bkgd_spectrum (:class:`~gdt.background.primitives.BackgroundSpectrum`):
                A modeled background spectrum
            distrib (str): The distribution from which the background is
                           simulated; either 'Poisson' or 'Gaussian'
        """
        if not isinstance(bkgd_spectrum, BackgroundSpectrum):
            raise TypeError('bkgd_spectrum must be a BackgroundSpectrum object')
        
        bkgd = BackgroundSpectrum(bkgd_spectrum.rates,
                                  bkgd_spectrum.rate_uncertainty,
                                  bkgd_spectrum.lo_edges,
                                  bkgd_spectrum.hi_edges,
                                  [self._sample_per] * bkgd_spectrum.size)

        if distrib == 'Poisson':
            self._spec_gen = VariablePoissonBackground(bkgd)
        elif distrib == 'Gaussian':
            self._spec_gen = VariableGaussianBackground(bkgd)
        else:
            raise ValueError("distrib can only be 'Poisson' or 'Gaussian'")

    def simulate(self, tstart, tstop):
        """Generate an EventList containing the individual counts from the 
        background simulation.
        
        Args:
            tstart (float): The start time of the simulation
            tstop (float): The stop time of the simulation
        
        Returns:
            (:class:`~gdt.core.data_primitives.EventList`)
        """
        # create the time grid
        dur = (tstop - tstart)
        numpts = int(round(dur / self._sample_per))
        time_array = np.linspace(tstart, tstop, numpts)

        # calculate the spectral amplitudes over the grid
        amps = self._time_func(time_array, *self._time_params)

        times = []
        chans = []
        for i in range(numpts):
            # update amplitude and generate the count spectrum
            self._spec_gen.amp = amps[i]
            self._event_gen.spectrum = self._whole_counts(
                next(self._spec_gen).counts)
            # generate the count arrival times for the time slice spectrum
            events = next(self._event_gen)
            if events is not None:
                times.extend((events[0] + time_array[i]).tolist())
                chans.extend(events[1].tolist())

        # create the eventlist
        ebounds = Ebounds.from_bounds(self._bkgd.lo_edges, self._bkgd.hi_edges)
        eventlist = EventList(times=times, channels=chans, ebounds=ebounds)
        eventlist.sort_time()
        return eventlist

    def to_tte(self, tstart, tstop, **kwargs):
        """Generate a TTE object containing the individual counts from the 
        background simulation.
        
        Args:
            tstart (float): The start time of the simulation
            tstop (float): The stop time of the simulation
            **kwargs: Options to pass to :class:`~gdt.core.tte.PhotonList`
        
        Returns:
            (:class:`~gdt.core.tte.PhotonList`)
        """
        eventlist = self.simulate(tstart, tstop)
        gti = Gti.from_bounds([eventlist.time_range[0]], 
                              [eventlist.time_range[1]])
        tte = PhotonList.from_data(eventlist, gti=gti, **kwargs)
        return tte

    def _whole_counts(self, counts):
        # because we can end up with fractional counts for the background
        # (the *rate* is what is typically modeled, and so no guarantee that
        #  counts will come out to whole integers)
        u = np.random.random(counts.size)
        whole_counts = counts.astype(int)
        mask = (counts - whole_counts) > u
        whole_counts[mask] += 1
        return whole_counts
