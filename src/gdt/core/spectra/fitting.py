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
import warnings
from datetime import datetime

import numpy as np
from numpy.random import multivariate_normal
from scipy.linalg import inv
from scipy.misc import derivative
from scipy.optimize import minimize, brentq
from scipy.stats import chi2

from gdt.core.background.primitives import BackgroundRates, BackgroundSpectrum
from gdt.core.data_primitives import EnergyBins
from gdt.core.pha import Bak

__all__ = ['SpectralFitter', 'SpectralFitterChisq', 'SpectralFitterCstat',
           'SpectralFitterPgstat', 'SpectralFitterPstat', 'chisq', 'cstat',
           'pgstat', 'pstat']


class SpectralFitter:
    """Base class for jointly-fitting spectra from multiple detectors. 
    This is a (large) wrapper around `scipy.optimize.minimize
    <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html>`_. 
    This class should not be instantiated by the user, but instead a child class 
    for a specific fit statistic should be instantiated. 
    
    Child classes should be designed with the following requirements:
    
    #. An ``__init__`` method that takes in the inputs in the Parent 
       ``__init__`` method. Can have additional arguments/keywords.
    #. A ``_eval_stat`` method that accepts the index into the data sets and
       the source model rate for that data set. It must return a single number 
       that represents the likelihood for a single detector/dataset.
    
    Parameters:
        pha_list (list of :class:`~gdt.core.pha.Pha`): 
            The PHA objects containg the count spectrum for each detector
        bkgd_list (list of :class:`~gdt.background.primitives.BackgroundRates`, \
                   list of :class:`~gdt.background.primitives.BackgroundSpectrum`, \
                   or list of :class:`~gdt.core.pha.Bak`): 
            The background rates object, background spectrum, or Bak object
            for each detector.  If given the background rates object, the times
            in the corresponding PHA object will be used for the limits of
            integration.
        rsp_list (list of :class:`~gdt.core.response.Rsp`): 
            The response object for each detector
        statistic (<function>): The function to be used as a fit statistic
        channel_masks (list of np.array(dtype=bool), optional): 
            A boolean mask for each detector that indicates which energy 
            channels are to be fit.
            
            Note:
                The channel_mask will typically be set by the 
                :attr:`~gdt.core.pha.Pha.valid_channels` attribute in the PHA 
                object.  You can set channel_masks to override the PHA 
                :attr:`~gdt.core.pha.Pha.valid_channels`, but this is usually 
                unnecessary.
        
        method (str, optional): 
            The fitting algorithm, which should be one of the options for 
            scipy.optimize.minimize. 
            
            Note:
                All solvers, with the exception of 'dogleg' and 'trust-exact',
                are supported at this time.          
    """

    def __init__(self, pha_list, bkgd_list, rsp_list, statistic,
                 channel_masks=None, method='SLSQP'):

        # check that we have the right numbers of data, backgrounds, responses,
        # and fit masks
        self._num_sets = len(pha_list)
        if (len(bkgd_list) != self._num_sets) or (
                len(rsp_list) != self._num_sets):
            raise ValueError('Number of datasets, backgrounds, and responses must be the same')

        # set the energy channel masks
        if channel_masks is not None:
            if len(channel_masks) != self._num_sets:
                raise ValueError('If channel_masks is set, must be equal to the number of datasets')

        # store response, exposure, and channel masks
        self._rsp = rsp_list
        self._exposure = np.array([pha.exposure for pha in pha_list])
        if channel_masks is None:
            channel_masks = []
            for pha in pha_list:
                chan_mask = np.zeros(pha.num_chans, dtype=bool)
                chan_mask[pha.valid_channels] = True
                channel_masks.append(chan_mask)
        self._chan_masks = channel_masks

        # extract observed counts and apply channel mask
        self._data = self._apply_masks([pha.data.counts for pha in pha_list])

        # extract background rates and variances and apply channel masks
        self._back_rates = []
        self._back_var = []
        for bkgd, pha in zip(bkgd_list, pha_list):
            if isinstance(bkgd, BackgroundRates):
                bkgd_spec = bkgd.integrate_time(*pha.time_range)
            elif isinstance(bkgd, BackgroundSpectrum):
                bkgd_spec = bkgd
            elif isinstance(bkgd, Bak):
                bkgd_spec = bkgd.data
            else:
                raise ValueError('Unknown Background object')
            self._back_rates.append(bkgd_spec.rates)
            self._back_var.append(bkgd_spec.rate_uncertainty ** 2)
        self._back_rates = self._apply_masks(self._back_rates)
        self._back_var = self._apply_masks(self._back_var)

        # fitter/function info
        self._stat = statistic
        self._method = method
        self._function = None
        self._res = None

    @property
    def covariance(self):
        """(np.array): The covariance matrix of the fitted parameters.  
            
        Note:
            The assumption of negligible higher-order terms in the Taylor 
            expansion may not be valid and therefore the covariance matrix 
            may not be valid."""
        if self._res is not None:
            cov = inv(-self.hessian)
            return cov

    @property
    def detectors(self):
        """(list of str): The detectors used in the fit"""
        return [rsp.detector for rsp in self._rsp]

    @property
    def dof(self):
        """(int): The number of degrees of freedom"""
        free_params = np.sum(self._function.free)
        datapts = np.sum([chan_mask.sum() for chan_mask in self._chan_masks])
        return datapts - free_params

    @property
    def energy_range(self):
        """(float, float): The full energy range of the fit"""
        emin = np.concatenate(
            self._apply_masks([rsp.ebounds.low_edges() for rsp in self._rsp]))
        emax = np.concatenate(
            self._apply_masks([rsp.ebounds.high_edges() for rsp in self._rsp]))
        return np.min(emin), np.max(emax)

    @property
    def function_components(self):
        """(list of str): The names of the individual model components, if 
        applicable"""
        if self._function is not None:
            try:
                return self._function.names
            except AttributeError:
                # We will gracefully return none if the function doesn't define names
                pass

    @property
    def function_name(self):
        """(str): The name of the model function being fit"""
        if self._function is not None:
            return self._function.name

    @property
    def hessian(self):
        """(np.array): The Hessian matrix of the fitted parameters, for which 
        the covariance matrix is -inverse(Hessian).  
            
        Note:
            The assumption of negligible higher-order terms in the Taylor 
            expansion may not be valid and therefore the Hessian matrix may 
            not be valid."""
        if self._res is not None:
            return self._hessian(self.parameters, self._function)

    @property
    def jacobian(self):
        """(np.array): The Jacobian vector of the likelihood as a function of 
        the parameters."""
        if self._res is not None:
            return self._res.jac

    @property
    def message(self):
        """(str): The message from the fitter on convergence/non-convergence"""
        if self._res is not None:
            return self._res.message

    @property
    def num_components(self):
        """(int): Number of function components"""
        if self._function is not None:
            return self._function.num_components

    @property
    def num_sets(self):
        """(int): Number of datasets/detectors being fit"""
        return self._num_sets

    @property
    def parameters(self):
        """(np.array): The parameter values at the maximum likelihood"""
        if self._res is not None:
            return self._res.x

    @property
    def statistic(self):
        """(float): The value of the fit statistic"""
        if self._res is not None:
            return self._res.fun

    @property
    def success(self):
        """(bool): True if the fit was successful, False otherwise."""
        if self._res is not None:
            return self._res.success

    @property
    def symmetric_errors(self):
        """(np.array): The 1-sigma uncertainty on the parameters assuming 
        higher-order terms of the Taylor expansion of the likelihood about the 
        maximum likelihood solution is negligible. If those terms are not 
        negligible, errors could be very incorrect or event NaN. See attributes 
        :attr:`covariance` and :attr:`hessian`."""
        if self._res is not None:
            id_array = np.identity(self.parameters.size, dtype=bool)
            return np.sqrt(np.abs(self.covariance[id_array]))

    def asymmetric_errors(self, cl=0.683):
        """Calculate the asymmetric uncertainty on each fitted parameter.
        This function uses `scipy.optimize.brentq
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html>`_
        to find the root of:
        
        crit = -2(LL_max-LL),
        
        where *crit* is the critical value defined by the number of free 
        parameters and the confidence level, and *LL_max* is the maximum 
        log-likelihood.  The brentq function uses a bracketing range to find the
        roots, and this function uses the best-fit values of the parameters for
        one side of the bracket and the min/max values defined in the fitted
        function for the other side of the bracket.
        
        Args:
            cl (float, optional): The confidence level of the uncertainties. 
                                  Default is 0.683.
        
        Returns:        
            (np.array): A (``nparams``, 2) array where each row reprensents the
                        (negative, positive) uncertainty on each parameter.
        """
        # must have a fit
        if self._res is None:
            raise RuntimeError('Fit has not been performed')

        if cl < 0.0 or cl > 1.0:
            raise ValueError('confidence level must be between 0-1')

        # calculate the critical value
        nparams = self.parameters.size
        crit = chi2.ppf(cl, nparams)

        # this is the objective function for which we are finding the root
        def the_func(param_v, index):
            params = self.parameters.copy()
            params[index] = param_v
            test_like = self._fold_model(self._function.fit_eval, params)
            delta = 2.0 * (-self._res.fun - test_like)

            return delta - crit

        # Find the lower and upper bound of the uncertainty.  Here, we use
        # Brent's Method to determine the bounds, which requires a bracketing
        # range of the solution.  One side of the bracket is the maximum 
        # likelihood value.  If the function parameter has a finite support,
        # the corresponding end of the support is used for the other side of the
        # bracket, otherwise, the other side of the bracket is iteratively
        # decreased (increased) until the lower (upper) bound is constrained
        # by Brent's Method.  This has to be done because we don't know a priori
        # how far away from maximum likelihood the bound is.

        # lower bound
        root1 = np.zeros(nparams)
        for i in range(nparams):
            param = self.parameters[i]
            while True:
                if self._function.min_values[i] == -np.inf:
                    if param < 0.0:
                        minval = 2.0 * param
                    elif param > 0.1:
                        minval = 0.5 * param
                    else:
                        minval = param - 1.0
                else:
                    minval = self._function.min_values[i]

                # if we're at the bounds, we can't bracket the peak
                try:
                    r, o = brentq(the_func, self.parameters[i], minval, args=(i,),
                                  full_output=True, maxiter=1000)
                    break
                except ValueError:
                    if minval == self._function.min_values[i]:
                        warnings.warn("Parameter exists at its lower bound")
                        r = minval
                        break
                param = minval

            root1[i] = r

        # upper bound
        root2 = np.zeros(nparams)
        for i in range(nparams):
            param = self.parameters[i]
            while True:
                if self._function.max_values[i] == np.inf:
                    if param < -0.1:
                        maxval = 0.5 * param
                    elif param > 0.0:
                        maxval = 2.0 * param
                    else:
                        maxval = 1.0 + param
                else:
                    maxval = self._function.max_values[i]

                try:
                    r, o = brentq(the_func, self.parameters[i], maxval, args=(i,),
                                  full_output=True, maxiter=1000)
                    break
                except ValueError:
                    # if we're at the bounds, we can't bracket the peak
                    if maxval == self._function.max_values[i]:
                        warnings.warn("Parameter exists at its upper bound")
                        r = maxval
                        break
                param = maxval

            root2[i] = r

        # convert the bounds to -/+ uncertainties
        errs = np.zeros((nparams, 2))
        errs[:, 0] = np.abs(self.parameters - root1)
        errs[:, 1] = np.abs(self.parameters - root2)

        return errs

    def data_count_spectrum(self, upper_limits_sigma=2.0):
        """The observed source count spectrum, which is (total - background).
        Data bins with source counts less than the model variance are converted
        to upper limits.
        
        Args:
            upper_limits_sigma (float, optional): The upper limits sigma. 
                                                  Default is 2. 
        
        Returns:
            (tuple): A 5-tuple containing, for each dataset:
            
            - *list of np.array*: Energy centroids
            - *list of np.array*: Energy channel half widths (in log-space)
            - *list of np.array*: differential count spectrum
            - *list of np.array*: 1-sigma count spectrum errors
            - *list of np.array*: Upper limits boolean mask
            
            The upper limits mask is a boolean mask where True indicates the 
            element is an upper limit.
        """

        if upper_limits_sigma <= 0.0:
            raise ValueError('upper_limits_sigma must be > 0.0')

        chanwidths = self._apply_masks(
            [rsp.drm.channel_widths for rsp in self._rsp])
        mvar = self.model_variance()

        # calculate the differential count spectrum and upper limit masks
        src_spectra = []
        ulmasks = []
        for i in range(self.num_sets):
            src_counts = (self._data[i] - self._back_rates[i] * self._exposure[i]) / (self._exposure[i] * chanwidths[i])

            ulmask = src_counts < upper_limits_sigma * np.sqrt(mvar[i])
            src_counts[ulmask] = upper_limits_sigma * np.sqrt(mvar[i])[ulmask]
            src_spectra.append(src_counts)
            ulmasks.append(ulmask)

        ecents = self._apply_masks([rsp.drm.channel_centroids for rsp in self._rsp])
        elo = self._apply_masks([np.array(rsp.drm.ebounds.low_edges()) for rsp in self._rsp])
        ehi = self._apply_masks([np.array(rsp.drm.ebounds.high_edges()) for rsp in self._rsp])
        ewidths = [np.array([ecents[i] - elo[i], ehi[i] - ecents[i]]) for i in range(len(ecents))]
        src_errs = [np.sqrt(var) for var in mvar]

        return ecents, ewidths, src_spectra, src_errs, ulmasks

    def fit(self, function, **kwargs):
        """Fit the given the function to the data.
        
        Args:
            function (:class:`~gdt.spectra.functions.Function`): 
                The function object to use
            **kwargs: Keyword arguments passed to the fitter
        """
        self._function = function

        # treat the fixed parameters
        init_params = np.array(function.default_values)[function.free].tolist()

        # do the fit, also returns the jacobian
        if self._method == 'dogleg' or self._method == 'trust-exact':
            raise RuntimeError('dogleg and trust-exact solvers not yet supported')

        if self._method == 'trust-ncg' or self._method == 'trust-krylov':
            hess = '3-point'
        else:
            hess = None

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = minimize(self._eval_stat_jac, init_params, args=(function,),
                           method=self._method, bounds=function.parameter_bounds(),
                           jac=True, hess=hess, **kwargs)
        self._res = res

    def model_count_spectrum(self):
        """The fitted model count spectrum.
        
        Returns:
            (list of :class:`~gdt.core.data_primitives.EnergyBins`)
        """
        if self._res is None:
            raise RuntimeError('Fit has not been performed')

        ebins = []
        for i in range(self.num_sets):
            rates = self._rsp[i].drm.fold_spectrum(self._function.fit_eval,
                                                   self.parameters)
            counts = rates * self._exposure[i]
            emin = np.array(self._rsp[i].ebounds.low_edges())
            emax = np.array(self._rsp[i].ebounds.high_edges())

            ebin = EnergyBins(counts[self._chan_masks[i]],
                              emin[self._chan_masks[i]],
                              emax[self._chan_masks[i]], self._exposure[i])
            ebins.append(ebin)

        return ebins

    def model_variance(self):
        """The differential source model variance, accounting for the data 
        variance and background model variance.
        
        Returns:
            (list of np.array)
        """
        if self._res is None:
            raise RuntimeError('Fit has not been performed')

        chanwidths = self._apply_masks(
            [rsp.drm.channel_widths for rsp in self._rsp])

        # model variance: need positive data counts, background rate and 
        # variance and exposure
        mvar = []
        for i in range(self.num_sets):
            rates = self._rsp[i].drm.fold_spectrum(self._function.fit_eval, self.parameters)
            model_rate = rates[self._chan_masks[i]]
            model_rate[model_rate < 0.0] = 0.0
            mvar.append(self._back_var[i] / (chanwidths[i]) ** 2
                        + (model_rate / chanwidths[i] + self._back_rates[i] / chanwidths[i])
                        / (np.abs(self._exposure[i]) * chanwidths[i]))

        return mvar

    def residuals(self, sigma=True):
        """The fit residuals for each dataset in units of differential 
        counts or model sigma
        
        Args:
            sigma (bool, optional): If True, return in units of model sigma, 
                otherwise return in units of differential counts. 
                Default is True.
        
        Returns:
            (tuple): A 4-tuple containing:
            
            - *list of np.array*: Energy centroids     
            - *list of np.array*: Energy channel half-widths (in log-space)
            - *list of np.array*: Fit residuals
            - *list of np.array*: Residual uncertainties
        """
        if self._res is None:
            raise RuntimeError('Fit has not been performed')

        chanwidths = self._apply_masks(
            [rsp.drm.channel_widths for rsp in self._rsp])
        # calculate the count spectrum model
        model = self.model_count_spectrum()

        # residuals are the differential model counts subtracted from the
        # differential source counts above background
        resid = []
        for i in range(self.num_sets):
            back_rates = self._back_rates[i] / chanwidths[i]
            rates = self._data[i] / (self._exposure[i] * chanwidths[i])
            resid.append((rates - back_rates) - model[i].rates_per_kev)

        # can calculate the residuals as a function of the model uncertainty
        model_var = self.model_variance()
        if sigma:
            resid = [resid[i] / np.sqrt(model_var[i]) for i in
                     range(self.num_sets)]
            resid_err = [np.ones_like(r) for r in resid]
        else:
            resid_err = [np.sqrt(var) for var in model_var]

        # the energy centroids and widths
        ecents = self._apply_masks([rsp.drm.channel_centroids for rsp in self._rsp])
        elo = self._apply_masks([np.array(rsp.drm.ebounds.low_edges()) for rsp in self._rsp])
        ehi = self._apply_masks([np.array(rsp.drm.ebounds.high_edges()) for rsp in self._rsp])
        ewidths = [np.array([ecents[i] - elo[i], ehi[i] - ecents[i]]) for i in range(len(ecents))]

        return ecents, ewidths, resid, resid_err

    def sample_flux(self, erange, num_samples=1000, num_points=1000, **kwargs):
        """Produce samples of the energy-integrated photon or energy flux.
        This uses the covariance matrix and assumes a multivariate normal 
        distribuion for the parameters.  Caveat Emptor.
        
        Args:
            erange (float, float): The energy range over which to calculate the
                                   flux
            num_samples (int): The number of spectra to sample. Default is 1000.
            num_points (int): The number of grid points in energy. Default is 
                              1000.
            **kwargs: Options to be passed to 
                      :meth:`~gbm.spectra.functions.Function.integrate`.
        
        Returns:
            (np.array)
        """
        if self._res is None:
            raise RuntimeError('Fit has not been performed')

        params = self.sample_parameters(size=num_samples)
        flux = [
            self._function.integrate(param_vec, erange, num_points=num_points,
                                     **kwargs) for param_vec in params]

        return flux

    def sample_parameters(self, **kwargs):
        """Produce samples of the fitted parameters.  This uses the covariance
        matrix and assumes a multivariate normal distribuion.  Caveat Emptor.
        
        Args:
            **kwargs: Options (like size) to be passed to `numpy.random.multivariate_normal \
            <https://numpy.org/doc/stable/reference/random/generated/numpy.random.multivariate_normal.html>`_
        
        Returns:
            (np.array)
        """
        if self._res is None:
            raise RuntimeError('Fit has not been performed')

        samples = multivariate_normal(self.parameters, self.covariance,
                                      **kwargs)
        return samples

    def sample_spectrum(self, which, num_samples=1000, num_points=1000,
                        components=False):
        r"""Produce samples of the model photon, energy, or :math:`\nu F_\nu` 
        spectrum. This uses the covariance matrix and assumes a multivariate 
        normal distribuion for the parameters. Caveat Emptor.
        
        Args:
            which (str): Which spectrum to return. Either 'photon', 'energy', 
                         or 'nufnu'
            num_samples (int, optional): The number of spectrum to sample.
                                         Default is 1000.
            num_points (int, optional): The number of grid points in energy. 
                                        Default is 1000.
            components (bool, optional): 
                If True, calculate the spectrum for each individual model 
                components (if there are multiple components). Default is False.
        
        Returns:
            (tuple): A 2-tuple containing:
            
            - *np.array*: The energy grid at which the spectrum is calculated
            - *np.array*: Array of shape (``num_samples``, ``num_points``) if 
              there is only one component or 
              (``num_samples``, ``num_components``, ``num_points``) if there 
              are multiple components.
        """
        if self._res is None:
            raise RuntimeError('Fit has not been performed')

        params = self.sample_parameters(size=num_samples)
        grid = np.logspace(*np.log10(self.energy_range), num_points)
        if self.num_components > 1:
            model = np.array([self._function.fit_eval(param, grid, components=components) for param in params])
        else:
            model = np.array([self._function.fit_eval(param, grid) for param in params])

        if which == 'energy':
            model *= grid[np.newaxis, :]
        elif which == 'nufnu':
            model *= grid[np.newaxis, :] ** 2
        else:
            pass
        return grid, model

    def save(self, filename):
        """Saves the fitter object to a compressed binary numpy file. The full
        state of the fitter is serialized and saved, and it can be restored 
        using the :meth:`load` method.  
        
        Args:
            filename (str): The complete filename to save the fitter object to
        """
        now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        np.savez_compressed(filename, time=now, obj=self)

    def spectrum(self, which, num_points=1000, components=False):
        r"""The model photon, energy, or :math:`\nu F_\nu` spectrum.
        
        Args:
            which (str): Which spectrum to return. Either 'photon', 'energy', 
                         or 'nufnu' 
            num_points (int, optional): The number of grid points in energy. 
                                        Default is 1000.
            components (bool, optional): If True, calculate the spectrum for 
                                         each individual model components (if 
                                         there are multiple components). 
                                         Default is False.
        
        Returns:
            (tuple): A 2-tuple containing:
            
            - *np.array*: The energy grid at which the spectrum is calculated 
            - *np.array* or *list of np.array*: model spectrum values, or 
              spectrum values for each component
        """
        if self._res is None:
            raise RuntimeError('Fit has not been performed')

        grid = np.logspace(*np.log10(self.energy_range), num_points)
        if self.num_components > 1:
            pmodel = self._function.fit_eval(self.parameters, grid,
                                             components=components)
        else:
            pmodel = self._function.fit_eval(self.parameters, grid)

            # multiply by energy or energy^2 if returning energy or nufnu spectrum
        if which == 'energy':
            pmodel *= grid
        elif which == 'nufnu':
            pmodel *= grid ** 2
        else:
            pass

        return grid, pmodel

    @classmethod
    def load(cls, filename):
        """Loads a fitter object that has been saved to disk in a compressed
        numpy file.  Thanks to the cool properties of the numpy format, this
        method will load the complete state of the fitter when it was saved.
        
        Args:
            filename (str): The complete filename of the saved file
        """
        loaded = np.load(filename, allow_pickle=True)
        obj = loaded['obj'].reshape(-1)[0]
        print('Fit loaded from {}'.format(loaded['time']))
        return obj

    def _apply_masks(self, a_list):
        """Apply the fit masks to a list of arrays.
        
        Args:
            a_list (list of np.array): The list of arrays
        
        Returns:
            (list of np.array)
        """
        return [np.asarray(one_list)[one_mask] for one_list, one_mask in zip(a_list, self._chan_masks)]

    def _eval_stat(self, set_num, src_model):
        """Evaluate the statistic for a single set. This must be defined by the
        derived class
        
        Args:
            set_num (int): The index number of the set
            src_model (np.array): The source model rates for the set
        
        Returns:
            (float)
        """
        raise NotImplementedError

    def _eval_stat_jac(self, params, function):
        """Evaluate the statistic (and Jacobian of the statistic as a function
        of the model parameters).
        
        Args:
            params (list): The parameters values
            function (class:`~.functions.Function`) The function object to use
        
        Returns:
            (float, np.array)
        """
        stat = self._fold_model(function.fit_eval, params)
        jac = self._jacobian(params, function)
        return -stat, -jac

    def _fold_model(self, function, params):
        """Folds the model throught the spectrum and calculates the fit 
        statistic
        
        Note:
            This is an empty function in the base class, and the inherited class
            must define this.
        
        Args:
            function (:class:`~.functions.Function`): The function object to use
            params (list): The parameters values
        
        Returns:
            (float)
        """
        stat = np.zeros(self.num_sets)
        for i in range(self.num_sets):
            # fold model through response and convert to raw model counts
            model = self._rsp[i].drm.fold_spectrum(function, params, channel_mask=self._chan_masks[i])

            # perform likelihood calculation for one dataset
            stat[i] = self._eval_stat(i, model)

        return stat.sum()

    def _hessian(self, params, function):
        """Calculate the Hessian of the fit statistic as a function of the 
        model parameters.  This is evalated numerically using finite differences.
        
        Args:
            params (list): The parameters values
            function (class:`~.functions.Function`) The function object to use
        
        Returns:
            (np.array)
        """
        num_params = len(params)
        hess = np.zeros((num_params, num_params))

        for i in range(num_params):
            for j in range(num_params):
                # Hessian is symmetric, so we only have to calculate either the
                # upper or lower triangle
                if i > j:
                    hess[i, j] = hess[j, i]
                    continue
                params_temp = params.copy()
                if params[i] != 0:
                    dx1 = np.abs((1.0 + function.delta_rel[i]) * params[i] - params[i])
                else:
                    dx1 = np.abs((1.0 + function.delta_abs[i]) * params[i] - params[i])
                if params[j] != 0:
                    dx2 = np.abs((1.0 + function.delta_rel[j]) * params[j] - params[j])
                else:
                    dx2 = np.abs((1.0 + function.delta_abs[j]) * params[j] - params[j])

                hess[i, j] = self._mixed_pderiv(function.fit_eval, params_temp,
                                                i, j, dx1, dx2)

        return hess

    def _jacobian(self, params, function):
        """Calculate the Jacobian of the fit statistic as a function of the 
        model parameters.  This is evalated numerically using finite differences.
        
        Args:
            params (list): The parameters values
            function (class:`~.functions.Function`) The function object to use
        
        Returns:
            (np.array)
        """
        num_params = len(params)
        grad = np.zeros(num_params)
        # for each parameter, do partial derivative using finite differences
        for i in range(num_params):
            params_temp = params.copy()
            dx = np.abs((1.0 + function.delta_rel[i]) * params[i] - params[i])
            grad[i] = self._pderiv(function.fit_eval, params_temp, index=i,
                                   dx=dx)
        return grad

    def _mixed_pderiv(self, function, params, index1, index2, dx1, dx2):
        """Mixed second partial derivative of function with respect to its
        parameters. This calls _pderiv.
        
        Args:
            function (class:`~.functions.Function`) The function object to use
            params (list): The parameters values
            index1 (int): The index into the parameter list, representing the 
                       parameter for which the first partial derivative will be 
                       calculated.
            index2 (int): The index into the parameter list, representing the 
                       parameter for which the second partial derivative will be 
                       calculated.
            dx1 (float): The step size for the first parameter
            dx2 (flot): The step size for the second parameter
        
        Returns:
            (float)
        """

        def wraps(x):
            params[index2] = x
            the_args = [function, params]
            return self._pderiv(*the_args, index=index1, dx=dx1)

        return derivative(wraps, params[index2], dx=dx2)

    def _pderiv(self, function, params, index=0, **kwargs):
        """Partial derivative of function with respect to one of its parameters
        This is a wrapper around scipy.misc.derivative to allow a partial
        derivative to be calculated        
        
        Args:
            function (class:`~.functions.Function`) The function object to use
            params (list): The parameters values
            index (int, optional): 
                The index into the parameter list, representing the parameter 
                for which the partial derivative will be calculated. Default is 0
            **kwargs: Options to pass to scipy.misc.derivative 
                      (such as the ``dx`` argument)
        
        Returns:
            (float)
        """

        def wraps(x):
            params[index] = x
            the_args = [function, params]
            return self._fold_model(*the_args)

        return derivative(wraps, params[index], **kwargs)


class SpectralFitterChisq(SpectralFitter):
    """Jointly-fit datasets using Chi-Squared Likelihood.

    Parameters:
        pha_list (list of :class:`~gdt.core.pha.Pha`): 
            The PHA objects containg the count spectrum for each detector
        bkgd_list (list of :class:`~gdt.background.primitives.BackgroundRates`, \
                   list of :class:`~gdt.background.primitives.BackgroundSpectrum`, \
                   or list of :class:`~gdt.core.pha.Bak`): 
            The background rates object, background spectrum, or Bak object
            for each detector.  If given the background rates object, the times
            in the corresponding PHA object will be used for the limits of
            integration.
        rsp_list (list of :class:`~gdt.core.response.Rsp`): 
            The response object for each detector
        channel_masks (list of np.array(dtype=bool), optional): 
            A boolean mask for each detector that indicates which energy 
            channels are to be fit.
            
            Note:
                The channel_mask will typically be set by the 
                :attr:`~gdt.core.pha.Pha.valid_channels` attribute in the PHA 
                object.  You can set channel_masks to override the PHA 
                :attr:`~gdt.core.pha.Pha.valid_channels`, but this is usually 
                unnecessary.
        
        method (str, optional): 
            The fitting algorithm, which should be one of the options for 
            scipy.optimize.minimize. 
            
            Note:
                All solvers, with the exception of 'dogleg' and 'trust-exact',
                are supported at this time.          
    """

    def __init__(self, pha_list, bkgd_list, rsp_list, **kwargs):
        super().__init__(pha_list, bkgd_list, rsp_list, chisq, **kwargs)

    def _eval_stat(self, set_num, src_model):
        # perform chisq calculation for one dataset
        return self._stat(self._data[set_num], self._back_rates[set_num],
                          self._back_var[set_num], src_model,
                          self._exposure[set_num])


class SpectralFitterPgstat(SpectralFitter):
    """Jointly-fit datasets using Profile Gaussian likelihood (PG-Stat).

    Parameters:
        pha_list (list of :class:`~gdt.core.pha.Pha`): 
            The PHA objects containg the count spectrum for each detector
        bkgd_list (list of :class:`~gdt.background.primitives.BackgroundRates`, \
                   list of :class:`~gdt.background.primitives.BackgroundSpectrum`, \
                   or list of :class:`~gdt.core.pha.Bak`): 
            The background rates object, background spectrum, or Bak object
            for each detector.  If given the background rates object, the times
            in the corresponding PHA object will be used for the limits of
            integration.
        rsp_list (list of :class:`~gdt.core.response.Rsp`): 
            The response object for each detector
        channel_masks (list of np.array(dtype=bool), optional): 
            A boolean mask for each detector that indicates which energy 
            channels are to be fit.
            
            Note:
                The channel_mask will typically be set by the 
                :attr:`~gdt.core.pha.Pha.valid_channels` attribute in the PHA 
                object.  You can set channel_masks to override the PHA 
                :attr:`~gdt.core.pha.Pha.valid_channels`, but this is usually 
                unnecessary.
        
        method (str, optional): 
            The fitting algorithm, which should be one of the options for 
            scipy.optimize.minimize. 
            
            Note:
                All solvers, with the exception of 'dogleg' and 'trust-exact',
                are supported at this time.          
    """

    def __init__(self, pha_list, bkgd_list, rsp_list, back_exp=None, **kwargs):
        super().__init__(pha_list, bkgd_list, rsp_list, pgstat, **kwargs)

        # background exposure (may == source exposure)
        if back_exp is None:
            back_exp = self._exposure
        self._back_exp = back_exp

    def _eval_stat(self, set_num, src_model):
        # perform pgstat calculation for one dataset
        return self._stat(self._data[set_num], src_model, self._exposure[set_num],
                          self._back_rates[set_num] * self._back_exp[set_num],
                          self._back_var[set_num], self._back_exp[set_num])


class SpectralFitterCstat(SpectralFitter):
    """Jointly-fit datasets using C-Stat (Poisson source with Poisson background)
    
    Parameters:
        pha_list (list of :class:`~gdt.core.pha.Pha`): 
            The PHA objects containg the count spectrum for each detector
        bkgd_list (list of :class:`~gdt.background.primitives.BackgroundRates`, \
                   list of :class:`~gdt.background.primitives.BackgroundSpectrum`, \
                   or list of :class:`~gdt.core.pha.Bak`): 
            The background rates object, background spectrum, or Bak object
            for each detector.  If given the background rates object, the times
            in the corresponding PHA object will be used for the limits of
            integration.
        rsp_list (list of :class:`~gdt.core.response.Rsp`): 
            The response object for each detector
        channel_masks (list of np.array(dtype=bool), optional): 
            A boolean mask for each detector that indicates which energy 
            channels are to be fit.
            
            Note:
                The channel_mask will typically be set by the 
                :attr:`~gdt.core.pha.Pha.valid_channels` attribute in the PHA 
                object.  You can set channel_masks to override the PHA 
                :attr:`~gdt.core.pha.Pha.valid_channels`, but this is usually 
                unnecessary.
        
        method (str, optional): 
            The fitting algorithm, which should be one of the options for 
            scipy.optimize.minimize. 
            
            Note:
                All solvers, with the exception of 'dogleg' and 'trust-exact',
                are supported at this time.          
    """

    def __init__(self, pha_list, bkgd_list, rsp_list, back_exp=None, **kwargs):
        super().__init__(pha_list, bkgd_list, rsp_list, cstat, **kwargs)

        # background exposure (may == source exposure)
        if back_exp is None:
            back_exp = self._exposure
        self._back_exp = back_exp

    def _eval_stat(self, set_num, src_model):
        # perform cstat calculation for one dataset
        return self._stat(self._data[set_num], src_model, self._exposure[set_num],
                          self._back_rates[set_num] * self._back_exp[set_num],
                          self._back_exp[set_num])


class SpectralFitterPstat(SpectralFitter):
    """Jointly-fit datasets using P-Stat (Poisson source with known background).
    This statistic assumes non-zero counts, therefore any channels with zero
    counts will be masked out and not used in fitting.

    Parameters:
        pha_list (list of :class:`~gdt.core.pha.Pha`): 
            The PHA objects containg the count spectrum for each detector
        bkgd_list (list of :class:`~gdt.background.primitives.BackgroundRates`, \
                   list of :class:`~gdt.background.primitives.BackgroundSpectrum`, \
                   or list of :class:`~gdt.core.pha.Bak`): 
            The background rates object, background spectrum, or Bak object
            for each detector.  If given the background rates object, the times
            in the corresponding PHA object will be used for the limits of
            integration.
        rsp_list (list of :class:`~gdt.core.response.Rsp`): 
            The response object for each detector
        channel_masks (list of np.array(dtype=bool), optional): 
            A boolean mask for each detector that indicates which energy 
            channels are to be fit.
            
            Note:
                The channel_mask will typically be set by the 
                :attr:`~gdt.core.pha.Pha.valid_channels` attribute in the PHA 
                object.  You can set channel_masks to override the PHA 
                :attr:`~gdt.core.pha.Pha.valid_channels`, but this is usually 
                unnecessary.
        
        method (str, optional): 
            The fitting algorithm, which should be one of the options for 
            scipy.optimize.minimize. 
            
            Note:
                All solvers, with the exception of 'dogleg' and 'trust-exact',
                are supported at this time.          
    """

    def __init__(self, pha_list, bkgd_list, rsp_list, **kwargs):
        super().__init__(pha_list, bkgd_list, rsp_list, pstat, **kwargs)

    def _eval_stat(self, set_num, src_model):
        # perform pstat calculation for one dataset
        return self._stat(self._data[set_num], src_model,
                          self._exposure[set_num], self._back_rates[set_num])


# --------------------------------------------------------------------------
# FIT STATISTICS

def chisq(obs_counts, back_rates, back_var, mod_rates, exposure):
    """Chi-Square Likelihood
    
    Args:
        obs_counts (np.array): The total observed counts
        back_rates (np.array): The background model count rates
        back_var (np.array): The background model count rate variance
        mod_rates (np.array): The model source rates
        exposure (float):  The source exposure
    
    Returns:    
        (float)
    """
    mask = (obs_counts - back_rates * exposure) > 0
    like = (((obs_counts[mask] - back_rates[mask] * exposure) - mod_rates[mask] * exposure) ** 2
            / (obs_counts[mask] + back_var[mask] * exposure ** 2))
    return -like.sum()


def cstat(obs_counts, mod_rates, exposure, back_counts, back_exp):
    """C-Statistic Likelihood for Poisson source + Poisson background.
    The "W" statistic from the `XSPEC Statistics Appendix
    <https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixStatistics.html>`_.
    
    Args:
        obs_counts (np.array): The total observed counts
        mod_rates (np.array): The model source rates
        exposure (float):  The source exposure
        back_counts (np.array): The background model counts
        back_exp (float): The background exposure
    
    Returns:
        (float)
    """
    sum_exp = exposure + back_exp
    # special cases for when the observed counts = 0 or when the background
    # model counts = 0
    mask1 = (obs_counts == 0)
    mask2 = (back_counts == 0) & (mod_rates < obs_counts / sum_exp)
    mask3 = (back_counts == 0) & (mod_rates >= obs_counts / sum_exp) & (mod_rates > 0)

    w = np.zeros_like(obs_counts)
    w[mask1] = exposure * mod_rates[mask1] - back_counts[mask1] * np.log(back_exp / sum_exp)
    w[mask2] = -back_exp * mod_rates[mask2] - obs_counts[mask2] * np.log(exposure / sum_exp)
    w[mask3] = (exposure * mod_rates[mask3] + obs_counts[mask3]
                * (np.log(obs_counts[mask3]) - np.log(exposure * mod_rates[mask3]) - 1.0))

    # for all other cases:
    mask = ~(mask1 | mask2 | mask3)
    mod_rates_nz = mod_rates[mask]
    obs_counts_nz = obs_counts[mask]
    back_counts_nz = back_counts[mask]

    d = np.sqrt(
        (sum_exp * mod_rates_nz - obs_counts_nz - back_counts_nz) ** 2 + 4.0 * sum_exp * back_counts_nz * mod_rates_nz)
    f = (obs_counts_nz + back_counts_nz - sum_exp * mod_rates_nz + d) / (2.0 * sum_exp)
    w[mask] = (exposure * mod_rates_nz + sum_exp * f - obs_counts_nz * np.log(exposure * mod_rates_nz + exposure * f)
               - back_counts_nz * np.log(back_exp * f) - obs_counts_nz * (1.0 - np.log(obs_counts_nz))
               - back_counts_nz * (1.0 - np.log(back_counts_nz)))

    return -w.sum()


def pstat(obs_counts, mod_rates, exposure, back_rates):
    """Likelihood for Poisson source + known background.
    The "pstat" statistic from the `XSPEC Statistics Appendix
    <https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixStatistics.html>`_.
    
    Note:
        Elements with zero counts are masked out and not figured in the 
        statistic.
    
    Args:
        obs_counts (np.array): The total observed counts
        mod_rates (np.array): The model source rates
        exposure (float):  The source exposure
        back_rates (np.array): The background model count rates
    
    Returns:    
        (float)
    """

    mask = (obs_counts > 0) & (mod_rates + back_rates > 0.0)
    pstat_val = (exposure * (mod_rates[mask] + back_rates[mask])
                 - obs_counts[mask] * np.log(exposure * (mod_rates[mask] + back_rates[mask]))
                 - obs_counts[mask] * (1.0 - np.log(obs_counts[mask])))
    return -pstat_val.sum()


def pgstat(obs_counts, mod_rates, exposure, back_counts, back_var, back_exp):
    """Profile Gaussian Likelihood. From the `XSPEC Statistics Appendix
    <https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixStatistics.html>`_.
    
    Args:
        obs_counts (np.array): The total observed counts
        mod_rates (np.array): The model source rates
        exposure (float): The source exposure
        back_counts (np.array): The background model counts
        back_var (np.array): The background model count variance
        back_exp (float): The background exposure
    
    Returns:    
        (float)
    """
    mask = (obs_counts > 0)
    pg = np.zeros_like(obs_counts)
    # special case for zero observed counts
    pg[~mask] = (exposure * mod_rates[~mask] + back_counts[~mask] * (exposure / back_exp)
                 - np.sqrt(back_var[~mask]) * (exposure / back_exp) ** 2 / 2.0)

    # for all other cases:
    obs_counts_nz = obs_counts[mask]
    mod_rates_nz = mod_rates[mask]
    back_counts_nz = back_counts[mask]
    back_var_nz = back_var[mask]

    d = np.sqrt((exposure * back_var_nz - back_exp * back_counts_nz + back_exp ** 2 * mod_rates_nz) ** 2
                - 4.0 * back_exp ** 2 * (exposure * back_var_nz * mod_rates_nz - obs_counts_nz * back_var_nz
                                         - back_exp * back_counts_nz * mod_rates_nz))

    f = ((-(exposure * back_var_nz - back_exp * back_counts_nz + back_exp ** 2 * mod_rates_nz) + d)
         / (2.0 * back_exp ** 2))

    pg[mask] = (back_exp * (mod_rates_nz + f) - obs_counts_nz * np.log(exposure * mod_rates_nz + exposure * f)
                + (back_counts_nz - back_exp * f) ** 2 / (2.0 * back_var_nz)
                - obs_counts_nz * (1.0 - np.log(obs_counts_nz)))

    return -pg.sum()
