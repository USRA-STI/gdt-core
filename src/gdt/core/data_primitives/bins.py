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
from .intervals import Ebounds

__all__ = ['Bins', 'ExposureBins', 'ChannelBins', 'TimeBins',
           'EnergyBins', 'TimeChannelBins', 'TimeEnergyBins']


class Bins():
    """A primitive class defining a set of histogram bins.
    The number of counts in each bin are assumed to be represented by a Poisson
    distribution unless the `count_uncerts` argument is specified with the 
    defined count uncertainties.
    
    Parameters:
        counts (np.array): The array of counts in each bin
        lo_edges (np.array): The low-value edges of the bins
        hi_edges (np.array): The high-value edges of the bins
        count_uncerts (np.array, optional): An array the same length as `counts`
                                            if the uncertainty is not Poisson.
    """
    def __init__(self, counts, lo_edges, hi_edges, count_uncerts=None):

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
        
        if count_uncerts is None:
            count_uncerts = np.sqrt(counts)
        try:
            iter(count_uncerts)
            self._count_uncerts = np.asarray(count_uncerts)
        except:
            raise TypeError('count_uncerts must be an iterable')
        if self._count_uncerts.size != self._counts.size:
            raise ValueError('count_uncerts must be same size as counts')
        
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
        return self._count_uncerts

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
    
    @classmethod
    def from_rates(cls, rates, rate_uncerts, lo_edges, hi_edges):
        """Create a Bins object from count rates and uncertainties.
        
        Args:
            rates (np.array): The count rates
            rate_uncerts (np.array): The count rate uncertainties
            lo_edges (np.array): The low-value edges of the bins
            hi_edges (np.array): The high-value edges of the bins
        
        Returns:
            (:class:`Bins`)
        """
        dt = hi_edges - lo_edges
        counts = rates * dt
        count_uncerts = rate_uncerts * dt
        return cls(counts, lo_edges, hi_edges, count_uncerts=count_uncerts)
    
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
        
        count_uncerts = self.count_uncertainty[mask]
        obj = Bins(self.counts[mask], self.lo_edges[mask], self.hi_edges[mask], 
                   count_uncerts=count_uncerts)
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
        count_uncerts (np.array, optional): An array the same length as `counts`
                                            if the uncertainty is not Poisson.
        precalc_good_segments (bool, optional): If True, calculates contiguous
                                                bin segments on initialization.
                                                Default is True.
    """
    def __init__(self, counts, lo_edges, hi_edges, exposure,
                 count_uncerts=None, precalc_good_segments=True):
        super().__init__(counts, lo_edges, hi_edges, count_uncerts=count_uncerts)

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

    @classmethod
    def from_rates(cls, rates, rate_uncerts, lo_edges, hi_edges, exposure, 
                   **kwargs):
        """Create an ExposureBins object from count rates and uncertainties.
        
        Args:
            rates (np.array): The count rates
            rate_uncerts (np.array): The count rate uncertainties
            lo_edges (np.array): The low-value edges of the bins
            hi_edges (np.array): The high-value edges of the bins
            exposure (np.array): The exposure of each bin
            precalc_good_segments (bool, optional): If True, calculates contiguous
                                                    bin segments on initialization.
                                                    Default is True.
        
        Returns:
            (:class:`ExposureBins`)
        """
        counts = rates * exposure
        count_uncerts = rate_uncerts * exposure
        return cls(counts, lo_edges, hi_edges, exposure, 
                   count_uncerts=count_uncerts, **kwargs)
    
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
        empty = self.__class__([], [], [], [], count_uncerts=[])

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
                new_counts, new_uncert, new_exposure, new_edges = method(histo.counts,
                                                                         histo.count_uncertainty,
                                                                         histo.exposure,
                                                                         edges, 
                                                                         *args)
                new_histo = self.__class__(new_counts, new_edges[:-1],
                                           new_edges[1:], new_exposure, 
                                           count_uncerts=new_uncert)
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
                             self.hi_edges[mask], self.exposure[mask], 
                             count_uncerts=self.count_uncertainty[mask])
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
        count_uncerts = histos[idx[0]].count_uncertainty
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
            count_uncerts = np.concatenate(
                (count_uncerts, histos[idx[i]].count_uncertainty[mask]))

        # new ExposureBins object
        merged_bins = cls(counts, lo_edges, hi_edges, exposure, 
                          count_uncerts=count_uncerts, **kwargs)
        return merged_bins

    @classmethod
    def sum(cls, histos):
        """Sum multiple ExposureBins together if they have the same bin edges.
        If the exposures are different between the histograms, they will be 
        averaged.
        
        Note:
            Count uncertainties are summed in quadrature.
        
        Args:
            histos (list of :class:`ExposureBins`):  
                A list containing the ExposureBins to be summed
        
        Returns:        
            (:class:`ExposureBins`)
        """
        counts = np.zeros(histos[0].size)
        count_var = np.zeros(histos[0].size)
        for histo in histos:
            assert histo.size == histos[0].size, \
                "The histograms must all have the same size"
            assert np.all(histo.lo_edges == histos[0].lo_edges), \
                "The histograms must all have the same support"
            counts += histo.counts
            count_var += histo.count_uncertainty ** 2
        
        count_uncerts = np.sqrt(count_var)
        
        # averaged exposure
        exposure = np.mean([histo.exposure for histo in histos], axis=0)

        sum_bins = cls(counts, histos[0].lo_edges, histos[0].hi_edges,
                       exposure, count_uncerts=count_uncerts)
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
    def create(cls, counts, chan_nums, exposure, continuous=True, **kwargs):
        """Create a :class:`ChannelBins` object from a list of channel numbers.
        
        Args:
            counts (np.array): The array of counts in each bin
            chan_nums (np.array): The energy channel numbers
            exposure (np.array): The exposure of each bin
            count_uncerts (np.array, optional): An array the same length as 
                                                `counts` if the uncertainty is 
                                                not Poisson.
            continuous (bool, optional): [Experimental] Whether the bins are continuous (meaning no gaps between channels).
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
        if continuous and len(chan_nums) > 1:
            hi_edges = np.concatenate((chan_nums[1:], [chan_nums[-1] + (chan_nums[1] - chan_nums[0])]))
        else:
            hi_edges = chan_nums + 1
        return cls(counts, chan_nums, hi_edges, exposure, **kwargs)

    @property
    def range(self):
        """(int, int): The channel range"""
        if self.size > 0:
            return (self.chan_nums[0], self.chan_nums[-1])

    @classmethod
    def from_rates(cls, rates, rate_uncerts, chan_nums, exposure, **kwargs):
        """Create an ChannelBins object from count rates and uncertainties.
        
        Args:
            rates (np.array): The count rates
            rate_uncerts (np.array): The count rate uncertainties
            chan_nums (np.array): The energy channel numbers
            exposure (np.array): The exposure of each bin
            precalc_good_segments (bool, optional): If True, calculates contiguous
                                                    bin segments on initialization.
                                                    Default is True.
        
        Returns:
            (:class:`ChannelBins`)
        """
        counts = rates * exposure
        count_uncerts = rate_uncerts * exposure
        return cls.create(counts, chan_nums, exposure, 
                          count_uncerts=count_uncerts, **kwargs)
    
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
                new_counts, new_uncert, new_exposure, new_edges = method(histo.counts,
                                                             histo.count_uncertainty,
                                                             histo.exposure,
                                                             edges, *args)
                new_histo = self.__class__(new_counts, new_edges[:-1],
                                           new_edges[1:], new_exposure, 
                                           count_uncerts=new_uncert)
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
                                    self.exposure[mask], 
                                    count_uncerts=self.count_uncertainty[mask])
        return obj


class TimeBins(ExposureBins):
    """A class defining a set of Time History (lightcurve) bins.

    Parameters:
        counts (np.array): The array of counts in each bin
        lo_edges (np.array): The low-value edges of the bins
        hi_edges (np.array): The high-value edges of the bins
        exposure (np.array): The exposure of each bin
        count_uncerts (np.array, optional): An array the same length as `counts`
                                            if the uncertainty is not Poisson.
        precalc_good_segments (bool, optional): If True, calculates contiguous
                                                bin segments on initialization.
                                                Default is True.
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
        count_uncerts (np.array, optional): An array the same length as `counts`
                                            if the uncertainty is not Poisson.
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

        Note:
            Count uncertainties are summed in quadrature.
                
        Args:
            histos (list of :class:`EnergyBins`):  
                A list containing the EnergyBins to be summed
        
        Returns:        
            (:class:`EnergyBins`)
        """
        counts = np.zeros(histos[0].size)
        uncerts = np.zeros(histos[0].size)
        exposure = 0.0
        for histo in histos:
            assert histo.size == histos[0].size, \
                "The histograms must all have the same size"
            assert np.all(histo.lo_edges == histos[0].lo_edges), \
                "The histograms must all have the same support"
            counts += histo.counts
            uncerts += histo.count_uncertainty ** 2
            exposure += histo.exposure

        uncerts = np.sqrt(uncerts)

        sum_bins = cls(counts, histos[0].lo_edges, histos[0].hi_edges,
                       exposure, count_uncerts=uncerts)
        return sum_bins


class TimeChannelBins():
    """A class defining a set of 2D Time/Energy-Channel bins.

    Parameters:
        counts (np.array): The array of counts in each bin
        tstart (np.array): The low-value edges of the time bins
        tstop (np.array): The high-value edges of the time bins
        exposure (np.array): The exposure of each bin
        chan_nums (np.array): The channel numbers in ascending order
        count_uncerts (np.array, optional): An array the same length as `counts`
                                            if the uncertainty is not Poisson.
        quality (np.array, optional): The spectrum quality flag
        precalc_good_segments (bool, optional): If True, calculates the 
                                                good time and channel segments 
                                                on initialization.    
    """
    def __init__(self, counts, tstart, tstop, exposure, chan_nums,
                 quality=None, count_uncerts=None, precalc_good_segments=True):

        try:
            iter(counts)
            self._counts = np.asarray(counts)
        except:
            raise TypeError('counts must be an iterable')
        if self._counts.ndim == 1:
            self._counts = self._counts.reshape(-1, 1)  
        elif self._counts.ndim > 2:
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

        if count_uncerts is None:
            count_uncerts = np.sqrt(self._counts)
        try:
            iter(count_uncerts)
            self._count_uncerts = np.asarray(count_uncerts)
        except:
            raise TypeError('count_uncerts must be an iterable')
        if self._count_uncerts.shape != self._counts.shape:
            raise ValueError('count_uncerts must be same shape as counts')

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
        return self._count_uncerts

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
                              ebounds.high_edges(), 
                              count_uncerts=self.count_uncertainty,
                              quality=self._quality)

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
        uncert = np.sqrt( (self.count_uncertainty[:, mask] ** 2).sum(axis=1) )

        obj = TimeBins(counts, self.tstart, self.tstop, self.exposure,
                       count_uncerts=uncert)
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
        counts = self.counts[mask, :].sum(axis=0)
        uncert = np.sqrt( (self.count_uncertainty[mask, :] ** 2).sum(axis=0) )
        exposure = np.full(counts.size, self.exposure[mask].sum())

        obj = ChannelBins.create(counts, self.chan_nums, exposure, 
                                 count_uncerts=uncert)
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
                new_uncert = []
                for i in range(num_times):
                    exposure = np.full(num_chans, histo.exposure[i])
                    new_cts, uncert, _, new_edges = method(histo.counts[i, :],
                                                   histo.count_uncertainty[i,:],
                                                   exposure, edges, *args)
                    new_counts.append(new_cts)
                    new_uncert.append(uncert)
                new_histo = TimeChannelBins(new_counts, bin.tstart, bin.tstop, 
                                            bin.exposure, new_edges[:-1], 
                                            count_uncerts=new_uncert)

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
                new_uncert = []
                for i in range(bin.num_chans):
                    new_cts, uncert, new_exposure, new_edges = method(
                        histo.counts[:, i],
                        histo.count_uncertainty[:, i],
                        histo.exposure,
                        edges, *args)
                    new_counts.append(new_cts)
                    new_uncert.append(uncert)
                new_counts = np.array(new_counts).T
                new_uncert = np.array(new_uncert).T

                new_histo = TimeChannelBins(new_counts, new_edges[:-1],
                                            new_edges[1:], new_exposure,
                                            bin.chan_nums, 
                                            count_uncerts=new_uncert)
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
                              count_uncerts=self.count_uncertainty[:, mask],
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
                  count_uncerts=self.count_uncertainty[mask, :],
                  quality=self.quality[mask])
        return obj

    @classmethod
    def from_rates(cls, rates, rate_uncerts, tstart, tstop, exposure, chan_nums,
                   **kwargs):
        """Create an TimeChannelBins object from count rates and uncertainties.
        
        Args:
            rates (np.array): The count rates
            rate_uncerts (np.array): The count rate uncertainties
            tstart (np.array): The low-value edges of the time bins
            tstop (np.array): The high-value edges of the time bins
            exposure (np.array): The exposure of each bin
            chan_nums (np.array): The channel numbers in ascending order
            quality (np.array, optional): The spectrum quality flag
            precalc_good_segments (bool, optional): If True, calculates the 
                                                    good time and energy 
                                                    segments on initialization.    
        
        Returns:
            (:class:`TimeChannelBins`)
        """
        counts = rates * exposure[:, np.newaxis]
        count_uncerts = rate_uncerts * exposure[:, np.newaxis]
        return cls(counts, tstart, tstop, exposure, chan_nums,  
                   count_uncerts=count_uncerts, **kwargs)
    
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
        uncert = histos[idx[0]].count_uncertainty
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
            uncert = np.hstack((uncert, histos[idx[i]].count_uncertainty[:, mask]))
            chan_nums = np.concatenate((chan_nums, histos[idx[i]].chan_nums[mask]))

        # new TimeChannelBins object
        merged_bins = cls(counts, tstart, tstop, exposure, chan_nums,
                          count_uncerts=uncert, quality=quality, **kwargs)
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
        uncert = histos[idx[0]].count_uncertainty
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
            uncert = np.vstack((uncert, histos[idx[i]].count_uncertainty[mask, :]))
            tstart = np.concatenate((tstart, histos[idx[i]].tstart[mask]))
            tstop = np.concatenate((tstop, histos[idx[i]].tstop[mask]))
            exposure = np.concatenate(
                (exposure, histos[idx[i]].exposure[mask]))
            quality = np.concatenate((quality, histos[idx[i]].quality[mask]))

        # new TimeChannelBins object
        merged_bins = cls(counts, tstart, tstop, exposure, chan_nums,
                          count_uncerts=uncert, quality=quality, **kwargs)
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
        count_uncerts (np.array, optional): An array the same length as `counts`
                                            if the uncertainty is not Poisson.
        quality (np.array, optional): The spectrum quality flag
        precalc_good_segments (bool, optional): If True, calculates the 
                                                good time and energy segments on
                                                initialization.    
    """
    def __init__(self, counts, tstart, tstop, exposure, emin, emax,
                 count_uncerts=None, quality=None, precalc_good_segments=True):

        try:
            iter(counts)
            self._counts = np.asarray(counts)
        except:
            raise TypeError('counts must be an iterable')
        if self._counts.ndim == 1:
            self._counts = self._counts.reshape(-1, 1)  
        elif self._counts.ndim > 2:
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

        if count_uncerts is None:
            count_uncerts = np.sqrt(self._counts)
        try:
            iter(count_uncerts)
            self._count_uncerts = np.asarray(count_uncerts)
        except:
            raise TypeError('count_uncerts must be an iterable')
        if self._count_uncerts.shape != self._counts.shape:
            raise ValueError('count_uncerts must be same shape as counts')

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
        return self._count_uncerts

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
        counts = self.counts[:, mask].sum(axis=1)
        uncert = np.sqrt( (self.count_uncertainty[:, mask] ** 2).sum(axis=1) )

        obj = TimeBins(counts, self.tstart, self.tstop, self.exposure,
                       count_uncerts=uncert)
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
        counts = self.counts[mask, :].sum(axis=0)
        uncert = np.sqrt( (self.count_uncertainty[mask, :] ** 2).sum(axis=0) )
        exposure = self.exposure[mask].sum()
        exposure = np.full(counts.size, exposure)

        obj = EnergyBins(counts, self.emin, self.emax, exposure, 
                         count_uncerts=uncert)
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
                new_uncert = []
                for i in range(num_times):
                    exposure = np.full(num_chans, histo.exposure[i])
                    new_cts, uncert, _, new_edges = method(histo.counts[i, :],
                                                   histo.count_uncertainty[i,:],
                                                   exposure, edges, *args)
                    new_counts.append(new_cts)
                    new_uncert.append(uncert)
                new_counts = np.array(new_counts)
                new_uncert = np.array(new_uncert)
                new_histo = TimeEnergyBins(new_counts, bin.tstart, bin.tstop,
                                           bin.exposure, new_edges[:-1],
                                           new_edges[1:],
                                           count_uncerts=new_uncert)
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
                new_uncert = []
                for i in range(bin.num_chans):
                    new_cts, uncert, new_exposure, new_edges = method(
                        histo.counts[:, i],
                        histo.count_uncertainty[:, i],
                        histo.exposure,
                        edges, *args)
                    new_counts.append(new_cts)
                    new_uncert.append(uncert)
                new_counts = np.array(new_counts).T
                new_uncert = np.array(new_uncert).T

                new_histo = TimeEnergyBins(new_counts, new_edges[:-1],
                                           new_edges[1:], new_exposure,
                                           bin.emin, bin.emax,
                                           count_uncerts=new_uncert)
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
                             count_uncerts=self.count_uncertainty[:, mask],
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
                  count_uncerts=self.count_uncertainty[mask, :],
                  quality=self.quality[mask])
        return obj

    @classmethod
    def from_rates(cls, rates, rate_uncerts, tstart, tstop, exposure, emin,
                   emax, **kwargs):
        """Create an TimeEnergyBins object from count rates and uncertainties.
        
        Args:
            rates (np.array): The count rates
            rate_uncerts (np.array): The count rate uncertainties
            tstart (np.array): The low-value edges of the time bins
            tstop (np.array): The high-value edges of the time bins
            exposure (np.array): The exposure of each bin
            emin (np.array): The low-value edges of the energy bins
            emax (np.array): The high-value edges of the energy bins
            quality (np.array, optional): The spectrum quality flag
            precalc_good_segments (bool, optional): If True, calculates the 
                                                    good time and energy 
                                                    segments on initialization.    
        
        Returns:
            (:class:`TimeEnergyBins`)
        """
        counts = rates * exposure[:, np.newaxis]
        count_uncerts = rate_uncerts * exposure[:, np.newaxis]
        return cls(counts, tstart, tstop, exposure, emin, emax,  
                   count_uncerts=count_uncerts, **kwargs)

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
        uncert = histos[idx[0]].count_uncertainty
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
            uncert = np.hstack((uncert, histos[idx[i]].count_uncertainty[:, mask]))
            emin = np.concatenate((emin, histos[idx[i]].emin[mask]))
            emax = np.concatenate((emax, histos[idx[i]].emax[mask]))

        # new TimeEnergyBins object
        merged_bins = cls(counts, tstart, tstop, exposure, emin, emax,
                          count_uncerts=uncert, quality=quality, **kwargs)
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
        uncert = histos[idx[0]].count_uncertainty
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
            uncert = np.vstack((uncert, histos[idx[i]].count_uncertainty[mask, :]))
            tstart = np.concatenate((tstart, histos[idx[i]].tstart[mask]))
            tstop = np.concatenate((tstop, histos[idx[i]].tstop[mask]))
            exposure = np.concatenate(
                (exposure, histos[idx[i]].exposure[mask]))
            quality = np.concatenate((quality, histos[idx[i]].quality[mask]))

        # new TimeEnergyBins object
        merged_bins = cls(counts, tstart, tstop, exposure, emin, emax,
                          count_uncerts=uncert, quality=quality, **kwargs)
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
