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
from gdt.core.data_primitives import TimeEnergyBins, EnergyBins, Gti, Ebounds

__all__ = ['BackgroundRates', 'BackgroundSpectrum']


class BackgroundRates(TimeEnergyBins):
    """Class containing the background rate data.

    Parameters:
        rates (np.array): The array of background rates in each bin
        rate_uncertainty (np.array): The array of background rate uncertainties 
                                     in each bin
        tstart (np.array): The low-value edges of the time bins
        tstop (np.array): The high-value edges of the time bins
        emin (np.array): The low-value edges of the energy bins
        emax (np.array): The high-value edges of the energy bins
        exposure (np.array, optional): The exposure of each bin    
    """
    def __init__(self, rates, rate_uncertainty, tstart, tstop, emin, emax,
                 exposure=None):
        
        try:
            iter(rates)
            rates = np.asarray(rates)
        except:
            raise TypeError('rates must be an iterable')
        if rates.ndim != 2:
            raise TypeError('rates must be a 2-dimensional array')

        try:
            iter(rate_uncertainty)
            rate_uncertainty = np.asarray(rate_uncertainty)
        except:
            raise TypeError('rate_uncertainty must be an iterable')
        if rate_uncertainty.ndim != 2:
            raise TypeError('rate_uncertainty must be a 2-dimensional array')
        if rate_uncertainty.shape != rates.shape:
            raise TypeError('rate_uncertainty must be the same shape as rates')
        
        if exposure is None:
            exposure = np.zeros_like(tstart)
        else:
            exposure = np.asarray(exposure)

        counts = np.squeeze(rates * exposure[:, np.newaxis])
        if counts.ndim == 1:
            counts = counts.reshape(tstart.size, emin.size)
        super().__init__(counts, tstart, tstop, exposure, emin, emax)
        self._count_uncertainty = np.squeeze(rate_uncertainty * exposure[:, np.newaxis])
        self._rates = rates.squeeze()
        self._rate_uncertainty = rate_uncertainty.squeeze()

    @property
    def count_uncertainty(self):
        """(np.array): The counts uncertainty in each bin"""
        return self._count_uncertainty

    @property
    def rate_uncertainty(self):
        """(np.array): The rate uncertainty in each bin"""
        return self._rate_uncertainty

    @property
    def rates(self):
        """(np.array): The rates in each Time-Energy Bin"""
        return self._rates

    def integrate_energy(self, emin=None, emax=None):
        """Integrate the over the energy axis.
        Limits on the integration smaller than the full range can be set.
        
        Args:
            emin (float, optional): The low end of the integration range. If not 
                               set, uses the lowest energy edge of the histogram
            emax (float, optional): The high end of the integration range. If not 
                              set, uses the highest energy edge of the histogram
        
        Returns:
            (:class:`BackgroundRates`)
        """
        if emin is None:
            emin = self.energy_range[0]
        if emax is None:
            emax = self.energy_range[1]

        mask = self._slice_energy_mask(emin, emax)
        emin = self.emin[mask][0]
        emax = self.emax[mask][-1]
        rates = np.nansum(self.rates[:, mask], axis=1).reshape(-1,1)
        rate_uncert = np.sqrt(
            np.nansum(self.rate_uncertainty[:, mask] ** 2, axis=1)).reshape(-1,1)

        obj = BackgroundRates(rates, rate_uncert, self.tstart, self.tstop,
                              np.array([emin]), np.array([emax]),
                              exposure=self.exposure)
        return obj

    def integrate_time(self, tstart=None, tstop=None):
        """Integrate the background over the time axis (producing a count rate
        spectrum). Limits on the integration smaller than the full range can 
        be set.
        
        Args:
            tstart (float, optional): The low end of the integration range.  
                          If not set, uses the lowest time edge of the histogram
            tstop (float, optional): The high end of the integration range. 
                         If not set, uses the highest time edge of the histogram
        
        Returns:
            (:class:`BackgroundSpectrum`)
        """
        if tstart is None:
            tstart = self.time_range[0]
        if tstop is None:
            tstop = self.time_range[1]

        mask = self._slice_time_mask(tstart, tstop)
        exposure = np.nansum(self.exposure[mask])
        rates = np.nansum(self.counts[mask, :], axis=0) / exposure
        rate_uncert = np.sqrt(np.nansum(self.count_uncertainty[mask, :] ** 2,
                                        axis=0)) / exposure
        exposure = np.full(rates.size, exposure)

        obj = BackgroundSpectrum(rates, rate_uncert, self.emin, self.emax,
                                 exposure)
        return obj

    def rebin_energy(self, method, *args, emin=None, emax=None):
        """Not implemented for BackgroundRates"""
        raise NotImplementedError

    def rebin_time(self, method, *args, tstart=None, tstop=None):
        """Not implemented for BackgroundRates"""
        raise NotImplementedError

    def slice_energy(self, emin, emax):
        """Perform a slice over an energy range and return a new BackgroundRates 
        object. Note that emin and emax values that fall inside a bin will 
        result in that bin being included.
        
        Args:
            emin (float): The start of the slice
            emax (float): The end of the slice
        
        Returns:           
            (:class:`BackgroundRates`)
        """
        mask = self._slice_energy_mask(emin, emax)
        cls = type(self)
        obj = cls(self.rates[:, mask], self.rate_uncertainty[:, mask], 
                  self.tstart, self.tstop, self.emin[mask], self.emax[mask], 
                  exposure=self.exposure)
        return obj

    def slice_time(self, tstart, tstop):
        """Perform a slice over a time range and return a new BackgroundRates 
        object. Note that tstart and tstop values that fall inside a bin will 
        result in that bin being included.
        
        Args:
            tstart (float): The start of the slice
            tstop (float): The end of the slice
        
        Returns:           
            (:class:`BackgroundRates`)
        """
        mask = self._slice_time_mask(tstart, tstop)
        cls = type(self)
        obj = cls(self.rates[mask, :], self.rate_uncertainty[mask,:], 
                  self.tstart[mask], self.tstop[mask], self.emin, self.emax, 
                  exposure=self.exposure[mask])
        return obj

    def to_bak(self, time_range=None, **kwargs):
        """Integrate over the time axis and produce a BAK object
        
        Args:
            time_range ((float, float), optional): 
                The time range to integrate over
            **kwargs: Options to pass to Bak.from_data()
        
        Returns:
            (:class:`~gdt.core.pha.Bak`)
        """
        from gdt.core.pha import Bak

        if time_range is None:
            time_range = self.time_range
        back_spec = self.integrate_time(*time_range)
        gti = Gti.from_list([time_range])
        
        bak = Bak.from_data(back_spec, gti=gti, **kwargs)
        return bak

    @classmethod
    def merge_time(cls, histos):
        """Merge multiple BackroundRates together along the time axis.
        
        Args:
            histos (list of :class:`BackgroundRates`): 
                A list containing the BackgroundRates to be merged
        
        Returns:
            (:class:`BackgroundRates`)
        """
        rates = np.vstack([histo.rates for histo in histos])
        rate_uncertainty = np.vstack([histo.rate_uncertainty \
                                      for histo in histos])
        bins = TimeEnergyBins.merge_time(histos)
        obj = cls(rates, rate_uncertainty, bins.tstart, bins.tstop,
                  bins.emin, bins.emax, exposure=bins.exposure)
        return obj

    @classmethod
    def sum_time(cls, bkgds):
        """Sum multiple BackgroundRates together if they have the same time 
        range.  Example use would be summing two backgrounds from two detectors.
        
        Args:
            bkgds (list of :class:`BackgroundRates`):
                A list containing the BackgroundRates to be summed
        
        Returns:
            (:class:`BackgroundRates`)
        """
        rates = np.zeros_like(bkgds[0].rates)
        rates_var = np.zeros_like(bkgds[0].rates)
        for bkgd in bkgds:
            assert bkgd.num_times == bkgds[0].num_times, \
                "The backgrounds must all have the same support"
            rates += bkgd.rates
            rates_var += bkgd.rate_uncertainty ** 2
            
        ebounds = Ebounds.from_bounds(bkgds[0].emin, bkgds[0].emax)
        for bkgd in bkgds[1:]:
            eb = Ebounds.from_bounds(bkgd.emin, bkgd.emax)
            ebounds = Ebounds.merge(ebounds, eb)

        # averaged exposure, sampling times
        exposure = np.mean([bkgd.exposure for bkgd in bkgds], axis=0)
        tstart = np.mean([bkgd.tstart for bkgd in bkgds], axis=0)
        tstop = np.mean([bkgd.tstop for bkgd in bkgds], axis=0)
        emin = ebounds.low_edges()
        emax = ebounds.high_edges()
        sum_bkgd = cls(rates, np.sqrt(rates_var), tstart, tstop, emin, emax,
                       exposure=exposure)
        return sum_bkgd


class BackgroundSpectrum(EnergyBins):
    """A class defining a Background Spectrum.

    Parameters:
        rates (np.array): The array of background rates in each bin
        rate_uncertainty (np.array): The array of background rate uncertainties 
                                     in each bin
        lo_edges (np.array): The low-value edges of the bins
        hi_edges (np.array): The high-value edges of the bins
        exposure (np.array): The exposure of each bin
    """
    def __init__(self, rates, rate_uncertainty, lo_edges, hi_edges, exposure):
        
        try:
            iter(rates)
            rates = np.asarray(rates)
        except:
            raise TypeError('rates must be an iterable')

        try:
            iter(rate_uncertainty)
            rate_uncertainty = np.asarray(rate_uncertainty)
        except:
            raise TypeError('rate_uncertainty must be an iterable')
        if rate_uncertainty.shape != rates.shape:
            raise TypeError('rate_uncertainty must be the same shape as rates')
        
        counts = rates * exposure
        super().__init__(counts, lo_edges, hi_edges, exposure)
        self._count_uncertainty = rate_uncertainty * exposure
        self._rates = rates
        self._rate_uncertainty = rate_uncertainty

    @property
    def count_uncertainty(self):
        """(np.array): The count uncertainty in each bin"""
        return self._count_uncertainty

    @property
    def rate_uncertainty(self):
        """(np.array): The count rate uncertainty of each bin"""
        return self._rate_uncertainty

    @property
    def rates(self):
        """(np.array): count rate of each bin"""
        return self._rates

    @classmethod
    def merge(cls, histos, **kwargs):
        """Not implemented for BackgroundSpectrum"""
        raise NotImplementedError
    
    def rebin(self, method, *args, emin=None, emax=None):
        """Not implemented for BackgroundSpectrum"""
        raise NotImplementedError

    def slice(self, emin, emax):
        """Perform a slice over an energy range and return a new 
        BackgroundSpectrum object. Note that the emin and emax values that fall 
        inside a bin will result in that bin being included.

        Args:
            emin (float): The low energy edge of the slice
            emax (float): The high energy of the slice
        
        Returns:
            (:class:`BackgroundSpectrum`)
        """
        emin_snap = self.closest_edge(emin, which='low')
        emax_snap = self.closest_edge(emax, which='high')

        mask = (self.lo_edges < emax_snap) & (self.hi_edges > emin_snap)
        obj = self.__class__(self.rates[mask], self.rate_uncertainty[mask],
                             self.lo_edges[mask], self.hi_edges[mask], 
                             self.exposure[mask])
        return obj

    @classmethod
    def sum(cls, histos):
        """Sum multiple BackgroundSpectrums together if they have the same 
        energy range (support).
        
        Args:
            histos (list of :class:`BackgroundSpectrum`):  
                A list containing the background spectra to be summed
        
        Returns:        
            (:class:`BackgroundSpectrum`)
        """
        counts = np.zeros(histos[0].size)
        count_variance = np.zeros(histos[0].size)
        exposure = 0.0
        for histo in histos:
            assert histo.size == histos[0].size, \
                "The histograms must all have the same size"
            assert np.all(histo.lo_edges == histos[0].lo_edges), \
                "The histograms must all have the same support"
            counts += histo.counts
            count_variance += histo.count_uncertainty**2
            exposure += histo.exposure
        
        rates = counts/exposure
        rate_uncertainty = np.sqrt(count_variance)/exposure
        
        sum_bins = cls(rates, rate_uncertainty, histos[0].lo_edges, 
                       histos[0].hi_edges, exposure)
        return sum_bins

