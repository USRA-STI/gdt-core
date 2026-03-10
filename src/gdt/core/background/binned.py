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
from scipy.interpolate import CubicSpline
from scipy.stats import norm
from scipy.optimize import curve_fit
from statsmodels.nonparametric.smoothers_lowess import lowess
from gdt.core.data_primitives import TimeChannelBins

__all__ = ['Polynomial', 'RoboLowess']


class Polynomial():
    """Class for performing a polynomial fit on Time-Energy data.
    The fit is performed over the time axis, treating each energy channel
    separately, although the fits are performed simultaneously.
    
    Parameters:
        counts (np.array): The array of counts in each bin, shape 
                           (``num_times``, ``num_chans``)
        tstart (np.array): The low-value edges of the time bins, shape 
                           (``num_times``,)
        tstop (np.array):  The high-value edges of the time bins, shape 
                           (``num_times``,)
        exposure (np.array): The exposure of each bin, shape (``num_times``,)
    """
    def __init__(self, counts, tstart, tstop, exposure):
        self._tstart = tstart
        self._tstop = tstop
        self._rate = counts / exposure[:, np.newaxis]
        self._livetime = exposure
        self._numtimes, self._numchans = self._rate.shape

        self._chisq = None
        self._dof = None
        self._order = None
        self._coeff = None
        self._covar = None

    @property
    def dof(self):
        """(np.array): The degrees-of-freedom for each channel"""
        return self._dof

    @property
    def statistic(self):
        """(np.array): The fit chi-squared statistic for each channel"""
        return self._chisq

    @property
    def statistic_name(self):
        """(str): 'chisq'"""
        return 'chisq'

    def fit(self, order=0):
        """Fit the data with a polynomial.
        Model variances are used for chi-squared via two fitting passes. 
        Adapted from the RMfit polynomial fitter.
        
        Args:
            order (int): The order of the polynomial
        
        Returns:        
            (np.array, np.array): The fitted model value and model uncertainty \
                                  at each input bin
        """
        assert order >= 0, 'Polynomial order must be non-negative'
        self._order = order
        self._coeff = np.zeros((order + 1, self._numchans))
        self._covar = np.zeros((order + 1, order + 1, self._numchans))

        # get basis functions and set up weights array
        tstart, tstop = self._tstart, self._tstop
        basis_func = self._eval_basis(tstart, tstop)
        weights = np.zeros((self._order + 1, self._numtimes, self._numchans))

        # Two-pass fitting
        # First pass uses the data variances calculated from data rates
        # Second pass uses the data variances calculated from model rates
        # 1) rate * livetime = counts
        # 2) variance of counts = counts
        # 3) variance of rate = variance of counts/livetime^2 = rate/livetime
        for iPass in range(2):
            if np.max(self._rate) <= 0.0:
                continue

            if iPass == 0:
                variance = self._rate / self._livetime[:, np.newaxis]
                idx = variance > 0.0
                for iCoeff in range(self._order + 1):
                    weights[iCoeff, idx] = basis_func[iCoeff, idx] / variance[
                        idx]
            else:
                variance = model / self._livetime[:, np.newaxis]
                idx = variance > 0.0
                if np.sum(idx) > 0:
                    for iCoeff in range(self._order + 1):
                        weights[iCoeff, idx] = basis_func[iCoeff, idx] / \
                                               variance[idx]
                else:
                    raise RuntimeError(
                        'SEVERE ERROR: Background model negative')

            # covariance matrix
            basis_func_list = basis_func.swapaxes(0,2).swapaxes(1,2) 
            weights_list = weights.swapaxes(0,2)
            rates_list = self._rate.T
            
            covar = np.array(
                list(map(np.dot, basis_func_list, weights_list))).T
            coeff = np.array(list(map(np.dot, rates_list, weights_list))).T

            if self._order >= 1:
                self._covar = np.linalg.inv(covar.T).T
            else:
                self._covar = 1.0 / covar

            # coefficients
            coeff_list = coeff.T            
            covar_list = self._covar.T
            coeff = np.array(list(map(np.dot, coeff_list, covar_list))).T

            # evaluate model
            self._coeff = coeff
            model = self._eval_model(tstart, tstop)
            

        # evaluate model uncertainty
        model_uncert = self._eval_uncertainty(tstart, tstop)

        # evaluate goodness-of-fit   
        self._chisq, self._dof = self._calc_chisq(model)

        return model, model_uncert

    def interpolate(self, tstart, tstop):
        """Interpolation of the fitted polynomial
        
        Args:
            tstart (np.array): The starting edge of each bin
            tstop (np.array): The ending edge of each bin
        
        Returns:        
            (np.array, np.array): The interpolated model value and model 
            uncertainty in each bin
        """
        interp = self._eval_model(tstart, tstop)
        interp_uncert = self._eval_uncertainty(tstart, tstop)
        return interp, interp_uncert

    def _calc_chisq(self, model):
        """Calculate the chi-squared goodness-of-fit for the fitted model.
        
        Args:
            model (np.array): The fitted model, shape (num_times, num_chans)

        Returns:
            (np.array, np.array) : The chi-squared goodness-of-fit and 
                                   degrees-of-freedom for each fitted channel
        """
        variance = model / self._livetime[:, np.newaxis]
        # do not calculate using bins with value <= 0.0
        idx = self._rate > 0.0
        chisq = [np.sum(
            (self._rate[idx[:, i], i] - model[idx[:, i], i]) ** 2 / variance[
                idx[:, i], i]) \
                 for i in range(self._numchans)]
        chisq = np.array(chisq)
        dof = np.sum(idx, axis=0) - (self._order + 1.0)
        return chisq, dof

    def _eval_basis(self, tstart, tstop):
        """Evaluates basis functions, which are the various polynomials 
        averaged over the time bins.
        
        Args: 
            tstart (np.array): The starting edge of each bin
            tstop (np.array): The ending edge of each bin
        
        Returns:        
            (np.array): The basis functions for each bin, shape \
                     (order + 1, num_times, num_chans)
        """
        numtimes = tstart.size
        dt = tstop - tstart
        zerowidth = (dt == 0.0)
        tstop[zerowidth] += 2e-6
        dt[zerowidth] = 2e-6

        basis_func = np.array(
            [(tstop ** (i + 1.0) - tstart ** (i + 1.0)) / ((i + 1.0) * dt) \
             for i in range(self._order + 1)])
        return np.tile(basis_func[:, :, np.newaxis], self._numchans)

    def _eval_model(self, tstart, tstop):
        """Evaluates the fitted model over the data
        
        Args: 
            tstart (np.array): The starting edge of each bin
            tstop (np.array): The ending edge of each bin
        
        Returns:        
            (np.array): The model value for each bin, shape \
                      (num_times, num_chans)
        """
        numtimes = tstart.size
        dt = tstop - tstart
        zerowidth = (dt == 0.0)
        tstop[zerowidth] += 2e-6
        dt[zerowidth] = 2e-6
        model = np.zeros((numtimes, self._numchans))
        for i in range(self._order + 1):
            model += (self._coeff[i, :] * (tstop[:, np.newaxis] ** (i + 1.0) - \
                       tstart[:, np.newaxis] ** (i + 1.0)) / \
                       ((i + 1.0) * dt[:, np.newaxis])).astype(float)
        
        if model.ndim == 1:
            model = model.reshape(-1, 1)
        return model

    def _eval_uncertainty(self, tstart, tstop):
        """Evaluates the uncertainty in the model-predicted values for the data 
        intervals based on the uncertainty in the model coefficients.
        
        Args: 
            tstart (np.array): The starting edge of each bin
            tstop (np.array): The ending edge of each bin
        
        Returns:        
            (np.array): The model uncertainty for each bin, shape \
                      (num_times, num_chans)
        """
        numtimes = tstart.size
        uncertainty = np.zeros((numtimes, self._numchans))
        basis_func = self._eval_basis(tstart, tstop)

        # formal propagation of uncertainty of fit coefficients to uncertainty of model
        covar_list = self._covar.T
        basis_func_list = basis_func.swapaxes(0,2).swapaxes(1,2)
        
        for i in range(numtimes):
            dot1 = np.array(list(map(np.dot, covar_list, 
                                     basis_func_list[:, :, i])))
            uncertainty[i, :] = np.array(list(map(np.dot, 
                                                  basis_func_list[:, :, i], dot1)))

        uncertainty = np.sqrt(uncertainty)
        return uncertainty

class RoboLowess:
    """Background fitting using LOWESS with iterative sigma-clipping.
    
    Args:
        counts (np.ndarray): Count array, shape (num_times, num_channels)
        tstart (np.ndarray): Time bin start edges, shape (num_times,)
        tstop (np.ndarray): Time bin stop edges, shape (num_times,)
        exposure (np.ndarray): Exposure of each bin, shape (num_times,)
        first_pass_chan_range (tuple, optional): (chan_min, chan_max) for Pass 1
    """

    def __init__(self, counts, tstart, tstop, exposure, 
                 *, first_pass_chan_range=None):
        
        num_chans = counts.shape[1]
        self._data = TimeChannelBins(counts, tstart, tstop, exposure, 
                                     np.arange(num_chans))
        
        if first_pass_chan_range is None:
            self._first_pass_chan_min = 0
            self._first_pass_chan_max = num_chans-1
        else:
            self._first_pass_chan_min = first_pass_chan_range[0]
            self._first_pass_chan_max = first_pass_chan_range[1]
        
        self._backgrounds = None
        self._pre_mask = None
        self._suppress_auto_refine = False
        self._refined_mask = None
        self._signal_mask = None
        self._refit_applied = False
        self._chisq = None
        self._dof = None
        self._statistic_name = 'chisq'

    @property
    def dof(self):
        """(np.array): The degrees-of-freedom for each channel"""
        return self._dof

    @property
    def statistic(self):
        """(np.array): The fit chi-squared statistic for each channel"""
        return self._chisq

    @property
    def statistic_name(self):
        """(str): 'chisq'"""
        return self._statistic_name


    def fit(self, win_size=None, temporal_resolution=None,
            min_win=0.4, max_win=0.95, spline_bc_type='clamped', lowess_iter=5):
        """Fit background model using two-pass LOWESS with sigma-clipping.
        
        Args:
            win_size (float, optional): LOWESS window fraction (0-1). 
                                       Auto-calculated if None.
            temporal_resolution (float, optional): Time resolution in seconds.
            min_win (float): Minimum window fraction (default 0.4)
            max_win (float): Maximum window fraction (default 0.95)
            spline_bc_type (str): Spline boundary condition (default 'clamped')
            lowess_iter (int): LOWESS robustness iterations (default 5)
        
        Returns:
            (np.ndarray, np.ndarray, list): Removed bin times, interpolated 
                                           background at removed times, 
                                           iteration diagnostics
        """
        # compute win_size
        if win_size is None:
            data_range = float(self._data.tstop[-1] - self._data.tstart[0])
            raw = 1.1 * (data_range ** -0.12) * (temporal_resolution ** -0.15) + 0.15
            win_size = min(max(min_win, raw), max_win)

        num_channels = self._data.rates.shape[1]
        backgrounds = np.zeros_like(self._data.rates)

        # Pass 1: Summed lightcurve to build background mask

        summed_data = self._data.integrate_channels(self._first_pass_chan_min, 
                                                    self._first_pass_chan_max)
        time_summed_all = np.asarray(summed_data.centroids).ravel()
        rate_summed_all = np.asarray(summed_data.rates).ravel()
        exp_summed = np.asarray(summed_data.exposure).ravel()

        if time_summed_all.size > 10:
            (times, bg_rates, bg_mask
             ) = self._residual_fit(
                 time_summed_all,
                 rate_summed_all,
                 exp_summed,
                 win_size=win_size,
                 lowess_iter=lowess_iter,
                 sigma_thresh=3.0,
                 return_mask=True,
                 use_pre_mask=True)

            # Background model on all summed times
            spline = CubicSpline(times, bg_rates,
                                 bc_type=spline_bc_type, extrapolate=True)
            bg_rates_full = spline(time_summed_all)
        else:
            # Fallback: everything is background with constant rate
            bg_mask = np.ones_like(time_summed_all, dtype=bool)
            bg_rates_full = np.full_like(rate_summed_all, 
                                         np.nanmean(rate_summed_all),
                                         dtype=float)

        # Track what we removed and interpolated background
        removed_bins_times = time_summed_all[~bg_mask]
        interp_back = (bg_rates_full[~bg_mask] *
                       exp_summed[~bg_mask])

        # Store for external access
        self._last_pass1_mask = bg_mask
        self._last_iteration_data = []

        # Background times from Pass 1
        bg_times_chan = time_summed_all[bg_mask]
        bg_exp_chan = exp_summed[bg_mask]

        # Pass 2: Per-channel refit using background times only

        for channel_index in range(num_channels):
            ch_rates = self._data.rates[bg_mask, channel_index]

            if bg_times_chan.size > 10:
                bg_times_ch, bg_rates_ch = self._residual_fit(
                                                bg_times_chan,
                                                ch_rates,
                                                bg_exp_chan,
                                                win_size=win_size,
                                                lowess_iter=lowess_iter,
                                                sigma_thresh=3.0,
                                                return_mask=False,
                                                use_pre_mask=False)

                # Spline to full time grid
                spline_channel = CubicSpline(
                                    bg_times_ch,
                                    bg_rates_ch,
                                    bc_type=spline_bc_type, 
                                    extrapolate=True)
                interp_bg_ch = spline_channel(
                                    self._data.time_centroids)
            else:
                # Not enough bins: use mean
                if ch_rates.size:
                    mean_channel_rate = np.nanmean(ch_rates)
                else:
                    mean_channel_rate = 0.0
                interp_bg_ch = np.full_like(
                    self._data.time_centroids, mean_channel_rate)

            backgrounds[:, channel_index] = interp_bg_ch

        # Store the full 2D background model (time × channel)
        self._backgrounds = backgrounds

        # Optional: Gaussian+Exponential refinement
        if not self._suppress_auto_refine:
            obs_counts = rate_summed_all * exp_summed
            bg_counts = bg_rates_full * exp_summed
            res_summed = (obs_counts - bg_counts) / np.sqrt(
                np.maximum(bg_counts, 1.0))

            # Check if refinement changes background
            bg_before = self._backgrounds.copy()

            refined_mask, signal_mask = self._refit_with_residual_model(
                                            residuals=res_summed, 
                                            threshold=0.5,
                                            win_size=win_size,
                                            spline_bc_type=spline_bc_type,
                                            lowess_iter=lowess_iter,
                                            first_pass_mask=bg_mask)
                                                                    
            if refined_mask is None:
                self._refit_applied = False
            else:
                self._refit_applied = not np.allclose(
                    self._backgrounds, bg_before)

        # Compute fit statistics on background bins (per channel)
        bg_mask = getattr(self, '_last_pass1_mask', None)
        if bg_mask is None:
            bg_mask = np.ones(self._data.time_centroids.size, dtype=bool)
        self._compute_statistics(bg_mask)

        return (removed_bins_times, interp_back,
                self._last_iteration_data)

    def _compute_statistics(self, bg_mask):
        num_times, num_channels = self._data.rates.shape
        counts = self._data.rates * self._data.exposure[:, None]
        expected = self._backgrounds * self._data.exposure[:, None]
        expected = np.maximum(expected, 1.0)

        counts_bg = counts[bg_mask]
        expected_bg = expected[bg_mask]

        resid = counts_bg - expected_bg
        self._chisq = np.sum((resid ** 2) / expected_bg, axis=0)
        dof_val = max(counts_bg.shape[0] - 1, 1)
        self._dof = np.full(num_channels, float(dof_val))


    def _residual_fit(self, time_all, rate_all, exposure_all,
                      win_size=0.5, lowess_iter=1, sigma_thresh=3.0,
                      core_snr_max=2.0, min_core_points=5,
                      return_mask=False, use_pre_mask=False, max_iter=10):
        """LOWESS with iterative sigma-clipping.
        
        Args:
            time_all (np.ndarray): Time values
            rate_all (np.ndarray): Count rates
            exposure_all (np.ndarray): Exposure values
            win_size (float): LOWESS window fraction (default 0.5)
            lowess_iter (int): LOWESS iterations (default 1)
            sigma_thresh (float): Threshold for outlier removal (default 3.0)
            core_snr_max (float): Max SNR for core distribution (default 2.0)
            min_core_points (int): Minimum points for core fit (default 5)
            return_mask (bool): Return background mask (default False)
            use_pre_mask (bool): Use pre-computed mask (default False)
            max_iter (int): Max iterations (default 10)
        
        Returns:
            If return_mask=False: (times, rates)
            If return_mask=True: (times, rates, mask)
        """
        time_all = np.asarray(time_all, dtype=float)
        rate_all = np.asarray(rate_all, dtype=float)
        exposure_all = np.asarray(exposure_all, dtype=float)
        num_bins = time_all.size

        # Not enough bins: flat background
        if num_bins < 10:
            mean_rate = np.nanmean(rate_all)
            flat = np.full_like(rate_all, mean_rate, dtype=float)
            if return_mask:
                return (time_all.copy(), flat,
                        np.ones_like(rate_all, dtype=bool))
            else:
                return time_all.copy(), flat

        if use_pre_mask and self._pre_mask is not None:
            background_mask = np.asarray(self._pre_mask, dtype=bool)
        else:
            background_mask = np.ones(num_bins, dtype=bool)

        for _ in range(max_iter):
            masked_times = time_all[background_mask]
            masked_rates = rate_all[background_mask]
            masked_exposure = exposure_all[background_mask]

            if masked_times.size < 10:
                break

            masked_times = np.asarray(masked_times).ravel()
            masked_rates = np.asarray(masked_rates).ravel()

            # LOWESS on rates
            lowess_result = lowess(masked_rates, masked_times,
                                   frac=win_size, it=lowess_iter)
            bg_rates_local = lowess_result[:, 1]

            # Residuals in counts space
            numerator = (masked_rates - bg_rates_local) * masked_exposure
            denominator = np.sqrt(
                np.maximum(bg_rates_local * masked_exposure, 1.0))
            snr = numerator / denominator

            # Core distribution
            core = snr[snr < core_snr_max]
            if core.size >= min_core_points:
                mu, sigma = norm.fit(core)
            else:
                mu = float(snr.mean())
                sigma_tmp = snr.std()
                sigma = float(sigma_tmp) if sigma_tmp > 0 else 1.0

            # Positive-sided clipping
            threshold = sigma_thresh * sigma
            keep_local = (snr < threshold)

            updated_mask = np.zeros_like(background_mask, dtype=bool)
            updated_mask[background_mask] = keep_local

            if np.array_equal(updated_mask, background_mask):
                background_mask = updated_mask
                break

            background_mask = updated_mask

        bg_times = time_all[background_mask]
        if bg_times.size < 2:
            mean_rate = np.nanmean(rate_all)
            if return_mask:
                return (time_all.copy(),
                        np.full_like(rate_all, mean_rate,
                                     dtype=float),
                        np.ones_like(rate_all, dtype=bool))
            else:
                return (time_all.copy(),
                        np.full_like(rate_all, mean_rate,
                                     dtype=float))

        bg_rates_input = rate_all[background_mask]
        bg_times = np.asarray(bg_times).ravel()
        bg_rates_input = np.asarray(bg_rates_input).ravel()

        # Final Lowess on the background subset
        lowess_result_final = lowess(bg_rates_input, bg_times,
                                     frac=win_size, it=lowess_iter)
        bg_rates = lowess_result_final[:, 1]

        if return_mask:
            return bg_times, bg_rates, background_mask
        else:
            return bg_times, bg_rates
    
    def _refit_with_residual_model(self, residuals, threshold=0.7,
                                   win_size=0.5, spline_bc_type='clamped',
                                   lowess_iter=1, first_pass_mask=None):
        """Refine background mask using Gaussian+Exponential model.
        
        Args:
            residuals (np.array): Residuals from Pass 1 fit
            threshold (float): Signal probability threshold (default 0.7)
            win_size (float): LOWESS window fraction (default 0.5)
            spline_bc_type (str): Spline BC (default 'clamped')
            lowess_iter (int): LOWESS iterations (default 1)
            first_pass_mask (np.array): Initial background mask
        
        Returns:
            (np.array, np.array): Refined mask, full signal mask
        """
        residuals = np.array(residuals)
        time_centroids = self._data.time_centroids

        # Identify finite residuals
        finite_mask = np.isfinite(residuals)

        # Combine with pass-1 mask
        if first_pass_mask is None:
            raise ValueError("Pass-1 mask required for refit")
        mask = finite_mask & first_pass_mask

        # Apply combined mask
        residuals = residuals[mask]
        time_centroids = time_centroids[mask]

        if residuals.size == 0:
            return None, None
        
        # Histogram of residuals
        num_bins = max(20, int(np.sqrt(len(residuals)))) 
        bins = np.linspace(residuals.min(), residuals.max(), num_bins + 1)
        hist, bin_edges = np.histogram(residuals, bins=bins, density=True)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        # Define Gaussian+Exponential model
        def gauss_exp(x, A, mu, sigma, B, lamb):
            gauss = A * norm.pdf(x, mu, sigma)
            exp_part = (B * lamb *
                        np.exp(np.clip(-lamb * (x - mu), -700, 700)) *
                        (x > mu))
            return gauss + exp_part

        try:
            popt, _ = curve_fit(gauss_exp, bin_centers, hist,
                                p0=[1.0, 0.0, 1.0, 0.1, 1.0])
        except RuntimeError:
            msg = ("Note: Gaussian+Exponential refinement failed. "
                   "Using Pass 2 background model.")
            print(msg)
            return None, None

        A, mu, sigma, B, lamb = popt

        # Compute signal probabilities
        gauss_pdf = A * norm.pdf(residuals, mu, sigma)
        exp_term = -lamb * (residuals - mu)
        exp_pdf = (B * lamb * np.exp(exp_term) * (residuals > mu))
        total_pdf = gauss_pdf + exp_pdf
        p_signal = np.divide(exp_pdf, 
                             total_pdf,
                             out=np.zeros_like(exp_pdf),
                             where=total_pdf > 0)

        # Define new mask: exclude high-probability signal bins
        signal_mask = p_signal > threshold

        # Map refined decisions back to full-length mask
        full_mask = first_pass_mask.copy()
        full_mask[mask] = ~signal_mask 

        # Refit using refined mask (Pass 1 and Pass 2 using new background mask)
        self._pre_mask = full_mask
        self._suppress_auto_refine = True

        try:
            self.fit(win_size=win_size, spline_bc_type=spline_bc_type,
                     lowess_iter=lowess_iter)
        finally:
            # Clear flags so future fits are normal
            self._suppress_auto_refine = False
            self._pre_mask = None

        # Full-length signal mask
        signal_mask_full = np.zeros_like(self._data.time_centroids,
                                         dtype=bool)
        signal_mask_full[mask] = signal_mask

        return full_mask, signal_mask_full