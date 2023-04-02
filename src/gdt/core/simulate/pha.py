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
from gdt.core.data_primitives import TimeEnergyBins, EnergyBins, Gti
from gdt.core.pha import Pha, Bak
from gdt.core.phaii import Phaii
from gdt.core.response import Rsp
from .generators import *

__all__ = ['PhaSimulator']

class PhaSimulator:
    """Simulate PHA data given a modeled background spectrum, detector response,
    source spectrum, and exposure.

    Parameters:
        rsp (:class:`~gdt.core.response.Rsp`): A detector response object
        function (:class:`~gdt.spectra.functions.Function`):
            A photon model function
        params (iterable): The parameters for the function
        exposure (float): The source exposure
        bkgd (:class:`~gdt.background.primitives.BackgroundSpectrum`):
            A modeled background spectrum
        bkgd_distrib (str): The distribution from which the background is
                            simulated; either 'Poisson' or 'Gaussian'
    """
    def __init__(self, rsp, function, params, exposure, bkgd, bkgd_distrib):
        self._rsp = rsp
        self._function = function
        self._params = params
        self._exposure = exposure
        self._src_gen = SourceSpectrumGenerator(rsp, function, params,
                                                exposure)
        self.set_background(bkgd, bkgd_distrib)

    def set_background(self, bkgd, bkgd_distrib):
        """Set/change the background model.
        
        Args:
            bkgd (:class:`~gdt.background.primitives.BackgroundSpectrum`):
                A modeled background spectrum
            bkgd_distrib (str): The distribution from which the background is
                                simulated; either 'Poisson' or 'Gaussian'
        """
        if not isinstance(bkgd, BackgroundSpectrum):
            raise TypeError('bkgd must be a BackgroundSpectrum object')
        bkgd_spectrum = BackgroundSpectrum(bkgd.rates, bkgd.rate_uncertainty,
                                           bkgd.lo_edges, bkgd.hi_edges,
                                           [self._exposure] * bkgd.size)
        if bkgd_distrib == 'Poisson':
            self._bkgd_gen = PoissonBackgroundGenerator(bkgd_spectrum)
        elif bkgd_distrib == 'Gaussian':
            self._bkgd_gen = GaussianBackgroundGenerator(bkgd_spectrum)
        else:
            raise ValueError(
                "bkgd_distrib can only be 'Poisson' or 'Gaussian'")

    def set_rsp(self, rsp):
        """Set/change the detector response.
        
        Args:
            rsp (:class:`~gdt.core.response.Rsp`): A detector response object        
        """
        if not isinstance(rsp, Rsp):
            raise TypeError('rsp must be a Rsp object')
        self._src_gen = SourceSpectrumGenerator(rsp, self._function,
                                                self._params,
                                                self._exposure)
        self._rsp = rsp

    def set_source(self, function, params, exposure):
        """Set/change the source spectrum.
        
        Args:
            function (:class:`~gdt.spectra.functions.Function`):  
                A photon model function        
            params (iterable): The parameters for the function
            exposure (float): The source exposure
        """
        self._src_gen = SourceSpectrumGenerator(self._rsp, function, params,
                                                exposure)
        self._exposure = exposure
        if self._bkgd_gen.__class__.__name__.startswith('Poisson'):
            distrib = 'Poisson'
        else:
            distrib = 'Gaussian'
        self.set_background(self._bkgd_gen._bkgd, distrib)        
        self._function = function
        self._params = params

    def simulate_background(self, num_sims):
        """Generate simulations of the modeled background spectrum.
        
        Args:
            num_sims (int): Number of simulations
        
        Returns:
            (list of :class:`~gdt.background.primitives.BackgroundSpectrum`)
        """
        return [next(self._bkgd_gen) for i in range(num_sims)]

    def simulate_source(self, num_sims):
        """Generate simulations of the source spectrum.
        
        Args:
            num_sims (int): Number of simulations
        
        Returns:
            (list of :class:`~gdt.core.data_primitives.EnergyBins`)
        """
        return [next(self._src_gen) for i in range(num_sims)]

    def simulate_sum(self, num_sims):
        """Generate simulations of the background + source spectrum.
        
        Args:
            num_sims (int): Number of simulations
        
        Returns:
            (list of :class:`~gdt.core.data_primitives.EnergyBins`)
        """
        summed = [None] * num_sims
        for i in range(num_sims):
            bkgd = next(self._bkgd_gen)
            src = next(self._src_gen)

            # since background model is formed from a rate, the background
            # "counts" won't be integers.  So we use the fractional part as
            # a probability to determine if we round up or truncate.
            bkgd_counts = bkgd.counts
            bkgd_counts[bkgd_counts < 0] = 0
            bkgd_counts_int = bkgd_counts.astype(int)
            bkgd_counts_frac = bkgd_counts - bkgd_counts_int
            extra_counts = (np.random.random(
                bkgd_counts_frac.size) > bkgd_counts_frac)
            bkgd_counts_int += extra_counts.astype(int)

            counts = bkgd_counts_int + src.counts
            summed[i] = EnergyBins(counts, src.lo_edges, src.hi_edges,
                                   src.exposure)
        return summed

    def to_bak(self, num_sims, tstart=None, tstop=None, **kwargs):
        """Produce BAK objects from simulations.
        
        Args:
            num_sims (int): Number of simulations
            tstart (float, optional): The start time. If not set, then is zero.
            tstop (float, optional): Then end time. If not set, then is the 
                                     exposure.
            **kwargs: Options passed to :class:`~gdt.core.pha.Bak`
        
        Returns:
            (list of :class:`~gdt.core.pha.Bak`)
        """
        if tstart is None:
            tstart = 0.0
        if tstop is None:
            tstop = tstart + self._exposure

        baks = self.simulate_background(num_sims)
        gti = Gti.from_bounds([tstart], [tstop])
        baks = [Bak.from_data(bak, gti=gti, **kwargs) for bak in baks]
        return baks

    def to_pha(self, num_sims, tstart=None, tstop=None, **kwargs):
        """Produce PHA objects of the background + source from simulations.
        
        Args:
            num_sims (int): Number of simulations
            tstart (float, optional): The start time. If not set, then is zero.
            tstop (float, optional): Then end time. If not set, then is the 
                                     exposure.
            **kwargs: Options passed to :class:`~gdt.core.pha.Pha`
        
        Returns:
            (list of :class:`~gdt.core.pha.Pha`)
        """
        if tstart is None:
            tstart = 0.0
        if tstop is None:
            tstop = tstart + self._exposure

        phas = self.simulate_sum(num_sims)
        gti = Gti.from_bounds([tstart], [tstop])
        phas = [Pha.from_data(pha, gti=gti, **kwargs) for pha in phas]
        return phas

    def to_phaii(self, num_sims, bin_width=None, **kwargs):
        """Produce a PHAII object by concatenating the simulations.
        
        Args:
            num_sims (int): Number of simulations
            bin_width (float, optional): The width of each time bin.  Must be
                >= the exposure.  If not set, the is the exposure.
            **kwargs: Options passed to :class:`~gdt.core.phaii.Phaii`
        
        Returns:
            (:class:`~gdt.core.phaii.Phaii`)
        """
        if bin_width is None:
            bin_width = self._exposure
        if bin_width < self._exposure:
            raise ValueError('bin_width cannot be less than exposure')

        phas = self.simulate_sum(num_sims)
        counts = np.vstack([pha.counts for pha in phas])
        edges = np.arange(num_sims + 1) * bin_width
        data = TimeEnergyBins(counts, edges[:-1], edges[1:],
                              [self._exposure] * num_sims, phas[0].lo_edges,
                              phas[0].hi_edges)
        
        gti = Gti.from_bounds([0.0], [bin_width*num_sims])
        phaii = Phaii.from_data(data, gti=gti, **kwargs)
        return phaii
