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
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid
from .bins import Bins
from .intervals import Ebounds


__all__ = ['ResponseMatrix']


class ResponseMatrix():
    """A class defining a 2D detector response matrix, with an (input) photon
    axis and (output) channel axis.

    Parameters:
        matrix (np.array): The 2D matrix, of shape 
                           (num_photon_bins, num_channels)
        emin (np.array): The low edges of the (input) photon bins
        emax (np.array): The high edges of the (input) photon bins
        chanlo (np.array): The low edges of the (output) energy channels
        chanhi (np.array): The high edges of the (output) energy channels    
    """
    def __init__(self, matrix, emin, emax, chanlo, chanhi):
        try:
            iter(matrix)
            self._matrix = np.asarray(matrix)
        except:
            raise TypeError('matrix must be an iterable')
        if self._matrix.ndim != 2:
            raise TypeError('matrix must be a 2-dimensional array')

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

        try:
            iter(chanlo)
            self._chanlo = np.asarray(chanlo).flatten()
        except:
            raise TypeError('chanlo must be an iterable')
        try:
            iter(chanhi)
            self._chanhi = np.asarray(chanhi).flatten()
        except:
            raise TypeError('chanhi must be an iterable')

        if (self._emin.size != self._emax.size):
            raise ValueError('emin and emax must be the same length')

        if (self._chanlo.size != self._chanhi.size):
            raise ValueError('chanlo and chanhi must be the same length')

        if (self._matrix.shape[0] != self._emin.size) or \
           (self._matrix.shape[1] != self._chanlo.size):
            raise ValueError('matrix axis 0 must have same length as emin and '\
                             'emax.  matrix axis 1 must have same length as '\
                             'chanlo and chanhi.')

    @property
    def channel_centroids(self):
        """(np.array): The geometric mean of the energy channels"""
        return np.sqrt(self._chanlo * self._chanhi)

    @property
    def channel_widths(self):
        """(np.array): The energy channel widths"""
        return self._chanhi - self._chanlo

    @property
    def ebounds(self):
        """(:class:`Ebounds`): The energy bounds"""
        return Ebounds.from_bounds(self._chanlo, self._chanhi)

    @property
    def matrix(self):
        """(np.array): The raw matrix"""
        return self._matrix

    @property
    def num_chans(self):
        """(int): The number of energy channels"""
        return self._matrix.shape[1]

    @property
    def num_ebins(self):
        """(int): The number of photon bins"""
        return self._matrix.shape[0]

    @property
    def photon_bins(self):
        """(:class:`Ebounds`): The photon bins"""
        return Ebounds.from_bounds(self._emin, self._emax)

    @property
    def photon_bin_centroids(self):
        """(np.array): The geometric mean of the photon bins"""
        return np.sqrt(self._emin * self._emax)

    @property
    def photon_bin_widths(self):
        """(np.array): The photon bin widths"""
        return self._emax - self._emin

    def channel_effective_area(self):
        """Returns the effective area as a function of recorded channel energy
        (integrated over incident photon bins).
        
        Returns:
            (:class:`Bins`)
        """
        bins = Bins(self._matrix.sum(axis=0), self._chanlo, self._chanhi)
        return bins

    def effective_area(self, energy, interp_kind='linear'):
        """Calculate the effective area at a given energy or an array of
        energies.  This function interpolates the DRM (in log10(cm^2)) to
        provide the effective area for any energy within the photon energy
        bounds of the DRM.
        
        Args:
            energy (float or np.array): The photon energy or energies at which
                                        to calculate the effective area.
            interp_kind (str, optional): The kind of interpolation to be 
                                         passed to scipy.interp1d.  Default is
                                         'linear'.
        
        Returns:
            (np.array)
        """
        try:
            energy = np.asarray(energy, dtype=float)
        except:
            raise TypeError('energy must be a positive float')

        if np.any(energy < self.photon_bins.range[0]) or \
           np.any(energy > self.photon_bins.range[1]):
            raise ValueError('energy must be within the photon energy bounds ' \
                             'of the DRM')

        effarea = self.photon_effective_area()
        diff_effarea = effarea.counts

        # create the DRM interpolator
        x = np.append(self._emin, self._emax[-1])
        y = np.append(diff_effarea, diff_effarea[-1])
        nz_mask = (y > 0)
        effarea_interpolator = interp1d(x[nz_mask], np.log10(y[nz_mask]),
                                        kind=interp_kind,
                                        fill_value='extrapolate')

        # do the interpolation
        effarea = 10.0**effarea_interpolator(energy)
        return effarea

    def fold_spectrum(self, function, params, channel_mask=None):
        """Fold a photon spectrum through a DRM to get a count spectrum
        
        Args: 
            function (<function>): 
                A photon spectrum function.  The function must accept a list of 
                function parameters as its first argument and an array of photon 
                energies as its second argument.  The function must return the 
                evaluation of the photon model in units of ph/s-cm^2-keV.
            params (list of float): A list of parameter values to be passed to
                                   the photon spectrum function
            channel_mask (np.array, optional): 
                A boolean mask where True indicates the channel is to be used 
                for folding and False indicates the channel is to not be used 
                for folding.  If omitted, all channels are used.
        
        Returns:        
            (np.array)
        """
        drm = self._matrix
        if channel_mask is not None:
            drm = drm[:, channel_mask]

        # evaluate photon model
        photon_model = function(params, self.photon_bin_centroids)

        # fold photon model through DRM
        counts = np.dot(drm.T, photon_model * self.photon_bin_widths)

        return counts

    def photon_effective_area(self):
        """Returns the effective area as a function of incident photon energy
        (integrated over recorded energy channels).
        
        Returns:        
            (:class:`Bins`)
        """
        bins = Bins(self._matrix.sum(axis=1), self._emin, self._emax)
        return bins

    def rebin(self, factor=None, edge_indices=None):
        """Rebins the channel energy axis of a DRM and returns a new response
        object.  Rebinning can only be used to downgrade the channel resolution
        and is constrained to the channel edges of the current DRM. 
        
        Rebinning can be performed by either defining an integral factor of the
        number of current energy channels to combine (e.g. factor=4 for 128
        energy channels would result in 32 energy channels), or by defining an
        index array into the current channel edges that will define the new set
        of edges.
        
        Args:
            factor (int, optional): The rebinning factor. Must set either this
                                    or `edge_indices`
            edge_indices (np.array, optional): The index array that represents
                                               which energy edges should remain
                                               in the rebinned DRM.
        
        Returns:
            (:class:`ResponseMatrix`)
        """
        if (factor is None) and (edge_indices is None):
            raise ValueError('Either factor or edge_indices must be set')

        elif factor is not None:
            try:
                factor = int(factor)
            except:
                raise TypeError('factor must be a positive integer')
            if factor < 1:
                raise ValueError('factor must be a positive integer')
            if (self.num_chans % factor) > 0:
                raise ValueError('factor must be factor of numchans: '\
                                 '{}'.format(self.num_chans))
            chanlo = self._chanlo.reshape(-1, factor)[:,0]
            chanhi = self._chanhi.reshape(-1, factor)[:,-1]

        elif edge_indices is not None:
            try:
                edge_indices = np.asarray(edge_indices)
            except:
                raise TypeError('new_edges must be a list, tuple or array')
            edge_indices = np.unique(edge_indices.astype(int))
            if (edge_indices[0] < 0) or (edge_indices[-1] > self.num_chans):
                raise ValueError('edge_indices outside valid range of ' \
                                 '0-{}'.format(self.num_chans))

            old_edges = np.append(self._chanlo, self._chanhi[-1])
            new_edges = old_edges[edge_indices]
            chanlo = new_edges[:-1]
            chanhi = new_edges[1:]

        # the new effective area matrix will just be summed over the channel
        # axis for the bins we are combining
        new_matrix = np.zeros((self.num_ebins, chanlo.size))
        sidx = (self._chanlo[:,np.newaxis] == chanlo[np.newaxis,:]).nonzero()[0]
        eidx = (self._chanhi[:,np.newaxis] == chanhi[np.newaxis,:]).nonzero()[0]
        for i in range(chanlo.size):
            new_matrix[:,i] = self._matrix[:,sidx[i]:eidx[i]+1].sum(axis=1)

        obj = type(self)(new_matrix, self._emin, self._emax, chanlo, chanhi)
        return obj

    def resample(self, num_photon_bins=None, photon_bin_edges=None,
                 num_interp_points=20, interp_kind='linear'):
        """Resamples the incident photon axis of a DRM and returns a new 
        ResponseMatrix object.  Resampling can be used to downgrade the photon 
        energy resolution, upgrade the resolution, and/or redefine the edges of 
        the incident photon bins.  By definition, the resampling can only be 
        performed within the photon energy range of the current object.
        
        The process for resampling is to create a high-resolution grid for each
        new photon bin, interpolate the differential effective area onto that
        grid (interpolation is performed on log10(effective area)), and then
        integrate the differential effective area over the grid points in each
        bin.  The final effective area is then scaled by the ratio of the 
        initial number of photon bins to the requested number of photon bins.
        
        Args:
            num_photon_bins (int, optional): The number of photon bins in the
                                             new DRM. The bin edges will be 
                                             generated logarithmically. Only set
                                             this or `photon_bin_edges`.
            photon_bin_edges (np.array, optional): The array of photon bin edges.
                                                   Only set this or 
                                                   `num_photon_bins`
            num_interp_points (int, optional): The number of interpolation points
                                               used to integrate over for each
                                               new photon bin. Default is 20.
            interp_kind (str, optional): The kind of interpolation to be 
                                         passed to scipy.interp1d.  Default is
                                         'linear'.
        
        Returns:
            (:class:`ResponseMatrix`)
        """
        if (num_photon_bins is None) and (photon_bin_edges is None):
            raise ValueError('Either num_photon_bins or photon_bin_edges must' \
                             ' be set')
        elif num_photon_bins is not None:
            try:
                num_photon_bins = int(num_photon_bins)
                assert num_photon_bins > 0
            except:
                raise TypeError('num_photon_bins must be a positive integer')
            photon_bin_edges = np.geomspace(self._emin[0], self._emax[-1],
                                            num_photon_bins+1)

        elif photon_bin_edges is not None:
            try:
                photon_bin_edges = np.asarray(photon_bin_edges)
            except:
                raise TypeError('photon_bin_edges must be an iterable')
            badmask = (photon_bin_edges < self._emin[0]) | \
                      (photon_bin_edges > self._emax[-1])
            if badmask.sum() > 0:
                raise ValueError('photon_bin_edges are beyond valid range of '\
                                 '{0:.2f}-{1:.2f}'.format(self._emin[0],
                                                          self._emax[-1]))
            photon_bin_edges = np.sort(photon_bin_edges)
            num_photon_bins = photon_bin_edges.size-1

        # differential effective area
        diff_effarea = self._matrix/self.photon_bin_widths[:,np.newaxis]

        # create an interpolation array over each new bin
        photon_bounds = np.vstack((photon_bin_edges[:-1], photon_bin_edges[1:]))
        interp_arr = np.geomspace(*photon_bounds, num_interp_points+2, axis=0)

        # create the DRM interpolator
        x = np.append(self._emin, self._emax[-1])
        y = np.vstack((diff_effarea, diff_effarea[-1,:]))
        y[y == 0] = 1e-10
        effarea_interpolator = interp1d(x, np.log10(y), axis=0,
                                        kind=interp_kind,
                                        fill_value='extrapolate')

        # interpolate the DRM over the photon interpolation array
        diff_effarea_interp = 10.**effarea_interpolator(interp_arr)
        badmask = np.isnan(diff_effarea_interp)
        diff_effarea_interp[badmask] = 0.0

        # integrate the DRM over the photon interpolation array
        diff_effarea_new = trapezoid(diff_effarea_interp,
                                 interp_arr[:,:,np.newaxis], axis=0)

        # scale by ratio of photon bins
        drm = diff_effarea_new / (self.num_ebins/num_photon_bins)

        # create response object
        obj = type(self)(drm, photon_bin_edges[:-1], photon_bin_edges[1:],
                         self._chanlo, self._chanhi)
        return obj

    def __repr__(self):
        return '<ResponseMatrix: {0} energy bins; {1} ' \
               'channels>'.format(self.num_ebins, self.num_chans)
