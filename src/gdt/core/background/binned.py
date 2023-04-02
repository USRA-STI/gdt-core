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

__all__ = ['Polynomial']


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

