import logging
logger = logging.getLogger(__name__)

from gdt.core.data_primitives import TimeRange, TimeBins
from gdt.core.background.binned import Polynomial
from gdt.core.binning.binned import rebin_by_time, rebin_by_edge_index
from gdt.misc.bayesian_blocks import bayesian_blocks

import matplotlib.pyplot as plt

from copy import copy

import numpy as np

from scipy.signal import peak_prominences, argrelmax

class BayesianBlocksLightcurve:
    """
    Get Bayesian blocks from a light curve while simultanously
    fitting the background and finding the signal time range.
    
    .. note::
       Call the method ``compute_bayesian_blocks()`` manually if you need
       to adjust any parameter of the algorithm. Otherwise it will be called with
       default parameters.
    
    Args:
        lc (TimeBins): Light curve
    """
    
    def __init__(self, lc):
        self._lc = lc

        self._lc_bayes = None
        self._lc_bb_index = None
        self._lc_bkg_counts = None
        self._bkg_times = None
        self._bkg_model = None    
        self._signal_range = None

        self._bb_args = None
        
    @property
    def lightcurve(self):
        """
        Get light curve

        Return:
            TimeBins
        """

        return self._lc

    def compute_bayesian_blocks(self,
                                p0 = 0.05, 
                                buffer_blocks = 5,
                                signal_range = None,
                                max_iter = 100):

        """
        Bayesian blocks with iterative background fitting
        
        Args:
            p0 (float): False alarm probability for bayesian blocks algorithm
            buffer_blocks (int): Define the exclusion zone around the signal
                to compute the background in multiple of the size of the first and last
                bayesian block that are part of the signal
            signal_range (TimeRange): Optional: provide a first guess for the time range
              that contained the signal.
           max_iter (int): Maximum number of iterations.

        Return:
           bb_index (array): Identified bayesian blocks, corresponding to indices of the ligh curve bins
           lc_bayes (TimeBins): Lightcurve bins in resulting bayesian blocks 

        """
        
        self._bb_args = {'p0':p0,
                         'buffer_blocks':buffer_blocks,
                         'signal_range':signal_range,
                         'max_iter':max_iter}
        
        # Standarize input
        lc = self.lightcurve
        
        if signal_range is None:
            bkg_times = [(lc.range[0], lc.range[-1])]
        else:
            if len(signal_range) != 2:
                raise ValueError("Wrong signal_range shape.")

            signal_range = np.sort(signal_range)

            bkg_times = [(lc.range[0], signal_range[0]),
                         (signal_range[-1], lc.range[-1])]

        # Iterative method
        previous_bkg_times = []
        
        for i in range(max_iter):

            # ---- Fit bkg ------
            lc_bkg_times = TimeBins.merge([lc.slice(ti,tf) for ti,tf in bkg_times])

            bkg_model = Polynomial(counts = lc_bkg_times.counts[:,np.newaxis], 
                                   tstart = lc_bkg_times.lo_edges, 
                                   tstop = lc_bkg_times.hi_edges, 
                                   exposure = lc_bkg_times.exposure)

            bkg_model.fit(order= min(2,i)) # Remove first mean and linear component

            bkg_rate, bkg_rate_err = bkg_model.interpolate(lc.lo_edges, lc.hi_edges)

            bkg_counts = bkg_rate.flatten() * lc.exposure

            # ---- Bayesian blocks --------

            # Make effective time bin ~bkg rate so the effective rate is
            # constant (homogeneous poisson process)
            # Same as Giacomo's "trick" in 3ML
            # https://github.com/threeML/threeML/blob/e31db70daf8777ce12be7aa694b21efc3f15dae0/threeML/utils/bayesian_blocks.py#L171  

            lc_eff = TimeBins(counts = lc.counts,
                              lo_edges = lc.lo_edges,
                              hi_edges = lc.hi_edges,
                              exposure = bkg_counts) 
            
            bb_index = bayesian_blocks(lc_eff, p0 = p0)

            lc_bayes = lc.rebin(rebin_by_edge_index, bb_index)

            lc_bkg_sub = TimeBins(lc.counts - bkg_counts, lc.lo_edges, lc.hi_edges, lc.exposure)
            lc_bayes_bkg_sub = lc_bkg_sub.rebin(rebin_by_edge_index, bb_index)

            # ----- Find peaks ----
            peaks = argrelmax(lc_bayes_bkg_sub.rates)[0]

            if len(peaks) == 0:
                # Need least 1 peak

                self._signal_range = signal_range
                self._bkg_times = bkg_times
                self._bkg_model = bkg_model
                self._lc_bkg_counts = bkg_counts
                self._lc_bb_index = bb_index
                self._lc_bayes = lc_bayes

                raise RuntimeError("Could not identify signal. Try with a different "
                                   "false positive rate or signal range initial guess.")

            prominence,left_base,right_base = peak_prominences(lc_bayes_bkg_sub.rates, peaks)

            # Start and end of signal based on peaks, in histogram bins
            leftmost_base = np.min(left_base)
            rightmost_base = np.max(right_base)

            new_start_signal = bb_index[leftmost_base + 1]
            new_stop_signal = bb_index[rightmost_base]

            signal_tstart = lc.lo_edges[new_start_signal]    
            signal_tstop  = lc.hi_edges[new_stop_signal-1]
            signal_range = (signal_tstart, signal_tstop)

            # Remove an extra chunk the size of the end blocks
            # Do not remove more than half the distance to the ends
            block_widths = lc_bayes.widths

            left_buffer = min(buffer_blocks*block_widths[left_base[0] + 1], (lc.lo_edges[new_start_signal] - lc.range[0])/2)
            right_buffer = min(buffer_blocks*block_widths[right_base[-1] - 1], (lc.range[1] - lc.hi_edges[new_stop_signal])/2)

            new_start,new_stop = np.digitize([lc.lo_edges[new_start_signal] - left_buffer,
                                              lc.hi_edges[new_stop_signal-1] + right_buffer], 
                                             lc.lo_edges) - 1

            # ----- Update background times ------
            # If they didn't change, then it has converge
            bkgex_tstart = lc.lo_edges[new_start]    
            bkgex_tstop  = lc.hi_edges[new_stop] 

            # Assume that the signal is fully contained. i.e. there is some background on both sides
            bkgex_tstart = max(bkgex_tstart, lc.lo_edges[bb_index[1]])
            bkgex_tstop = min(bkgex_tstop, lc.hi_edges[bb_index[-2]-1])
            
            bkg_times = [(lc.lo_edges[0], bkgex_tstart), (bkgex_tstop, lc.hi_edges[-1])]

            # Check if we are repeating a pattern.
            # Usually, when the method converges, 2 consecutive iteration have the same bayesian
            # blocks, but sometimes there is a 2-3 pattern that repeats with almost identical
            # block representations
            if bkg_times in previous_bkg_times and i >= 2:
                break

            previous_bkg_times += [bkg_times]

            if i == max_iter - 1:
                logger.warn("Maximum number of iterations reached without converging.")
                
        # Cache
        self._signal_range = signal_range
        self._bkg_times = bkg_times
        self._bkg_model = bkg_model
        self._lc_bkg_counts = bkg_counts
        self._lc_bb_index = bb_index
        self._lc_bayes = lc_bayes
        
    def _compute_bayesian_blocks(self):
        """
        Call compute_bayesian_blocks with default parameters if it hasn't been called before
        """
        if self._lc_bayes is None:
            self.compute_bayesian_blocks()            
    
    def plot(self, ax = None, rebin_dt = None, subtract_bkg = False):
        """
        Quick plot of original light curve, Bayesian blocks, estimated background, and
        identified signal range.

        Parameters:
        ax (matplotlib.axes): Axes where to draw the plot. A new 
                one will be created by default.
        rebin_dt (float): Rebin the raw light curve in bins of this width
        
        """
        self._compute_bayesian_blocks()
        
        if ax is None:
            fig,ax = plt.subplots()

        # Underying data
        lc = self.lightcurve

        if subtract_bkg:
            lc_bkg_sub = TimeBins(lc.counts - self.bkg_counts, lc.lo_edges, lc.hi_edges, lc.exposure)

            if rebin_dt is not None:
                lc_bkg_sub = lc_bkg_sub.rebin(rebin_by_time, rebin_dt)

        if rebin_dt is not None:
            lc_rebin = lc.rebin(rebin_by_time, rebin_dt)
        else:
            lc_rebin = lc

        if subtract_bkg:
            lc_plotting = lc_bkg_sub
        else:
            lc_plotting = lc_rebin

        # Plots
        ax.errorbar(lc_plotting.centroids, lc_plotting.rates,
                    xerr = [lc_plotting.centroids - lc_plotting.lo_edges,
                            lc_plotting.hi_edges - lc_plotting.centroids],
                    yerr = lc_rebin.rate_uncertainty,
                    ls = 'none', color = '.7',
                    label = 'Raw data')

        # Bayes on top
        lc_bayes = self.bb_lightcurve

        if subtract_bkg:
            bkg_bb = np.diff(np.append(0, np.cumsum(self.bkg_counts))[self._lc_bb_index])
            bb_bkg_sub_counts = self.bb_lightcurve.counts - bkg_bb
            bb_eff_rate = bb_bkg_sub_counts / self.bb_lightcurve.exposure

            ax.plot(np.append(lc_bayes.lo_edges, lc_bayes.hi_edges[-1]),
                    np.append(bb_eff_rate, bb_eff_rate[-1]),
                    drawstyle='steps-post',
                    label='Bayesian blocks')

        else:

            # Weight them by the background, since they are changes on top of the
            # non-homogeneous background

            bb_eff_rate = self.bkg_counts / lc.exposure

            for n, (start_idx, stop_idx) in enumerate(zip(self._lc_bb_index[:-1], self._lc_bb_index[1:])):
                # Make bb_eff_rate have the same average the actual data by adding a constant rate
                bb_eff_rate[start_idx:stop_idx] += lc_bayes.rates[n] - np.mean(bb_eff_rate[start_idx:stop_idx])

            ax.plot(np.append(lc.lo_edges, lc.hi_edges[-1]),
                    np.append(bb_eff_rate, bb_eff_rate[-1]),
                    drawstyle='steps-post',
                    label='Bayesian blocks')


        # Background
        if not subtract_bkg:
            ax.plot(lc.centroids, self.bkg_counts/lc.exposure, color = 'red', ls = ':',
                    label = "Fitted background")

        # Vertical lines showing the start and stop of the identified signal
        ax.axvline(self.signal_range.tstart, ls = "--", color = 'olive', label = "Signal start/stop")
        ax.axvline(self.signal_range.tstop, ls = "--", color = 'olive')

        ax.legend()

        ax.set_xlabel("Time")

        if subtract_bkg:
            ax.set_ylabel("Background-subtracted rate [Hz]")
        else:
            ax.set_ylabel("Rate [Hz]")

        return ax

    @property
    def bb_lightcurve(self):
        """
        Get the light curve rebin in bayesian blocks

        Return:
            
        """

        self._compute_bayesian_blocks()
        
        return self._lc_bayes
        
    @property
    def bkg_counts(self):
        """
        Return the estimated background counts corresponding to each bin in the
        original light curve.

        Return:
           array
        """

        self._compute_bayesian_blocks()
        
        return self._lc_bkg_counts

    @property
    def bkg_model(self):
        """
        Fitted background model

        Return:
           Binned background model. Currently only Polynomial 
        """
        return self._bkg_model
    
    @property
    def bkg_times(self):
        """
        Time ranges used to fit the background

        Return:
            list of TimeRange
        """
        return [TimeRange(*tr) for tr in self._bkg_times]

    @property
    def signal_range(self):
        """
        Get the time range of identified signal

        Return:
           TimeRange
        """

        self._compute_bayesian_blocks()

        return TimeRange(*self._signal_range)

    def duration(self, quantile):
        """
        Compute the symmetric time range that contains a given quantile of bkg-subtracted counts

        Args:
            quantile (float or array): Quantile(s) e.g. .9 for T90, [1,.9,.5] for [T100, T90, T50]

        Return:
           float or array: Same shape as quantile
        """

        self._compute_bayesian_blocks()

        lc = self.lightcurve
        bkg_counts = self.bkg_counts
        signal_range = self.signal_range

        # Percentile calculation
        start_signal,stop_signal = np.digitize([signal_range.tstart, signal_range.tstop],
                                               lc.lo_edges)-1

        cumtime = lc.hi_edges[start_signal:stop_signal] - lc.lo_edges[start_signal]

        # Background subtraction 
        cumcounts = np.cumsum(lc.counts[start_signal:stop_signal] -
                              bkg_counts[start_signal:stop_signal])


        half_inv_quant = (1-np.array(quantile))/2

        tquant = np.diff(np.interp([half_inv_quant, 1-half_inv_quant],
                                    cumcounts/cumcounts[-1],
                                    cumtime), axis = 0)
        
        if np.isscalar(quantile):
            tquant = tquant.item()
        else:
            tquant = tquant.reshape(np.shape(quantile))

        return tquant

    def duration_error(self, quantile, containment = .68, nsamples = 100):
        """
        Compute the confidence interval for the burst duration using random sampling.

        Args:
            quantile (float or array): Quantile(s) e.g. .9 for T90, [1,.9,.5] for [T100, T90, T50]
            containment (float or array): Confidence interval containment fraction

        Return:
            array: Negative and positive error. Shaped (2, size(containment), size(quantile)) 
        """
        
        self._compute_bayesian_blocks()

        lc = self.lightcurve
        lc_bayes = self._lc_bayes
        bkg_counts = self._lc_bkg_counts
        signal_range = self._signal_range

        # Create fake LC with mean signal counts on top of calculated background
        bb_mean_counts = lc_bayes.rates[np.digitize(lc.centroids, lc_bayes.lo_edges)-1] * lc.exposure
        
        start_signal,stop_signal = np.digitize(signal_range, lc.lo_edges)-1

        mean_counts = copy(bkg_counts)
        mean_counts[start_signal:stop_signal] = bb_mean_counts[start_signal:stop_signal]
        
        lc_sim = copy(lc)
        
        # Fluctuate counts and get samples
        tquant_samples = np.empty((nsamples,) + np.shape(quantile))

        for s in range(nsamples):

            logger.debug(f"Sample {s}/{nsamples}")

            lc_sim = TimeBins(counts = np.random.poisson(mean_counts),
                              lo_edges = lc_sim.lo_edges,
                              hi_edges = lc_sim.hi_edges,
                              exposure = lc_sim.exposure) 

            bb_lc_sim = self.__class__(lc_sim)

            bb_lc_sim.compute_bayesian_blocks(**self._bb_args)

            tquant = bb_lc_sim.duration(quantile)
                        
            tquant_samples[s] = tquant

        # Quantiles, subtracting the median
        containment = np.asarray(containment)

        tquant_median = np.quantile(tquant_samples,
                                    q = .5,
                                    axis = 0)
        
        tquant_bounds = np.quantile(tquant_samples,
                                    q = [(1-containment)/2, 1-(1-containment)/2],
                                    axis = 0)

        tquant_error = tquant_bounds - tquant_median
        
        # Add the size of the bins to the error
        binning_error = (lc.widths[start_signal] + lc.widths[stop_signal])/2
        tquant_error[0] -= binning_error
        tquant_error[-1] += binning_error
        
        return tquant_error

        
