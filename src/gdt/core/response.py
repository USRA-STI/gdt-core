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
import os
import astropy.io.fits as fits
import numpy as np
from scipy.interpolate import interp1d

from .data_primitives import ResponseMatrix, EnergyBins, TimeBins
from .file import FitsFileContextManager
from .headers import FileHeaders

__all__ = ['Rsp', 'Rsp2']

class Rsp(FitsFileContextManager):
    """Class for single-DRM responses.
    """
    def __init__(self):
        super().__init__()
        self._ebounds = None
        self._drm = None
        self._trigtime = None
        self._tstart = None
        self._tstop = None
        self._detector = None
        
        # used for compressing/decompressing matrix
        self._ngrp = None
        self._fchan = None
        self._nchan = None
    
    @property
    def detector(self):
        """(str): The name of the detector this response is associated with"""
        return self._detector
    
    @property
    def drm(self):
        """(:class:`~.data_primitives.ResponseMatrix`): The response matrix"""
        return self._drm
    
    @property
    def ebounds(self):
        """(:class:`~.data_primitives.Ebounds`): The energy channel bounds"""
        return self._ebounds
    
    @property
    def num_chans(self):
        """(int): Number of energy channels"""
        return self.drm.num_chans

    @property
    def num_ebins(self):
        """(int): Number of photon energy bins"""
        return self.drm.num_ebins

    @property
    def tcent(self):
        """(float): The center time of the interval over which the DRMs is 
        valid"""
        return (self.tstart + self.tstop) / 2.0

    @property
    def trigtime(self):
        """(float): The trigger time, if applicable"""
        return self._trigtime

    @property
    def tstart(self):
        """(float): The start time of the interval over which the DRM is valid"""
        return self._tstart

    @property
    def tstop(self):
        """(float): The end time of the interval over which the DRM is valid"""
        return self._tstop

    def fold_spectrum(self, function, params, channel_mask=None, exposure=1.0):
        """Fold a photon spectrum through a DRM to get a count spectrum.
        This function differs from 
        :meth:`~.data_primitives.ResponseMatrix.fold_spectrum` by
        taking in an optional exposure and returning an 
        :class:`~.data_primitives.EnergyBins` containing the count spectrum.
        
        Args: 
            function (<function>): 
                A photon spectrum function.  The function must accept a list of 
                function parameters as its first argument and an array of photon 
                energies as its second argument.  The function must return the 
                evaluation of the photon model in units of ph/s-cm^2-keV.
            params (list of float): A list of parameter values to be passed to
                                   the photon spectrum function
            channel_mask (np.array, optional): 
                A Boolean mask where True indicates the channel is to be used 
                for folding and False indicates the channel is to not be used 
                for folding.  If omitted, all channels are used.
            exposure (float, optional): The exposure in seconds. Default is 1.
        
        Returns:        
            (:class:`~.data_primitives.EnergyBins`)
        """
        try:
            exposure = float(exposure)
        except:
            raise TypeError('exposure must be a positive float')
        if exposure <= 0.0:
            raise ValueError('exposure must be positive')
        
        count_spec = self.drm.fold_spectrum(function, params, 
                                            channel_mask=channel_mask)
        counts = count_spec * exposure
        lo_edges = self.ebounds.low_edges()
        hi_edges = self.ebounds.high_edges()
        if channel_mask is not None:
            lo_edges = np.array(lo_edges)[channel_mask]
            hi_edges = np.array(hi_edges)[channel_mask]

        ebins = EnergyBins(counts, lo_edges, hi_edges, exposure)
        return ebins

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
                                    or ``edge_indices``.
            edge_indices (np.array, optional): The index array that represents
                                               which energy edges should remain
                                               in the rebinned DRM.
        
        Returns:
            (:class:`Rsp`)
        """
        new_drm = self.drm.rebin(factor=factor, edge_indices=edge_indices)
        tstart = self.tstart
        tstop = self.tstop
        if (tstart is not None) and (self.trigtime is not None):
            tstart += self.trigtime
        if (tstop is not None) and (self.trigtime is not None):
            tstop += self.trigtime
            
        headers = self._build_headers(new_drm.num_chans, new_drm.num_ebins)
        
        rsp = self.from_data(new_drm, start_time=tstart, stop_time=tstop, 
                             trigger_time=self.trigtime, headers=headers)
        rsp._ngrp = self._ngrp
        rsp._fchan = self._fchan
        rsp._nchan = self._nchan
        return rsp
    
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
                                             this or ``photon_bin_edges``.
            photon_bin_edges (np.array, optional): The array of photon bin edges.
                                                   Only set this or 
                                                   ``num_photon_bins``.
            num_interp_points (int, optional): The number of interpolation points
                                               used to integrate over for each
                                               new photon bin. Default is 20.
            interp_kind (str, optional): The kind of interpolation to be 
                                         passed to scipy.interp1d.  Default is
                                         'linear'.
        
        Returns:
            (:class:`Rsp`)
        """
        new_drm = self.drm.resample(num_photon_bins=num_photon_bins, 
                                    photon_bin_edges=photon_bin_edges,
                                    num_interp_points=num_interp_points,
                                    interp_kind=interp_kind)
        tstart = self.tstart
        tstop = self.tstop
        if (tstart is not None) and (self.trigtime is not None):
            tstart += self.trigtime
        if (tstop is not None) and (self.trigtime is not None):
            tstop += self.trigtime

        headers = self._build_headers(new_drm.num_chans, new_drm.num_ebins)
        
        rsp = self.from_data(new_drm, start_time=tstart, stop_time=tstop, 
                             trigger_time=self.trigtime, headers=headers)
        rsp._ngrp = self._ngrp
        rsp._fchan = self._fchan
        rsp._nchan = self._nchan
        return rsp
    
    @classmethod
    def from_data(cls, drm, filename=None, start_time=None, stop_time=None,
                  trigger_time=None, headers=None, detector=None):
        """Create a Rsp object from a 
        :class:`~.data_primitives.ResponseMatrix` object.
        
        Args:
            drm (:class:`~.data_primitives.ResponseMatrix`): The DRM
            filename (str, optional): The filename of the object
            start_time (float, optional): The start time for the response
            stop_time (float, optional): The stop time for the response
            trigger_time (float, optional): The trigger time, if applicable.
            headers (:class:`~.headers.FileHeaders`): The file headers
            detector (str, optional): The detector name
                 
        Returns:
            (:class:`Rsp`)
        """
        obj = cls()
        obj._filename = filename
        obj._detector = detector
        
        # set data and ebounds
        if not isinstance(drm, ResponseMatrix):
            raise TypeError('data must be of type ResponseMatrix')
        obj._drm = drm
        obj._ebounds = drm.ebounds
        
        # set time info
        if trigger_time is not None:
            if trigger_time < 0.0:
                raise ValueError('trigger_time must be non-negative')
        obj._trigtime = trigger_time
        obj._tstart = start_time
        if obj._trigtime is not None and obj._tstart is not None:
            obj._tstart -= obj._trigtime     
        obj._tstop = stop_time
        if obj._trigtime is not None and obj._tstop is not None:
            obj._tstop -= obj._trigtime     
                                
        # set headers
        if headers is not None:
            if not isinstance(headers, FileHeaders):
                raise TypeError('headers must be of type FileHeaders')
        obj._headers = headers        
        
        return obj

    def _build_headers(self, num_chans, num_ebins):
        """This builds the headers for the FITS file.  This method needs
        to be specified in the inherited class.  The method should construct
        the headers from the minimum required arguments and additional keywords
        and return a :class:`~.headers.FileHeaders` object.
        
        Args:
            num_chans (int): Number of detector energy channels
            num_ebins (int): Number of photon energy bins
        
        Returns:
            (:class:`~.headers.FileHeaders`)
        """
        pass

    def __repr__(self):
        s = '<{0}: '.format(self.__class__.__name__)
        if self.filename is not None:
            s += '{};'.format(self.filename)
        if self.trigtime is not None:
            s += '\n trigger time: {};'.format(self.trigtime)
        s += '\n time range ({0}, {1});'.format(self.tstart, self.tstop)
        s += '\n {0} energy bins; {1} channels>'.format(self.num_ebins, 
                                                        self.num_chans)
        return s
    

class Rsp2(FitsFileContextManager):
    """Class for multiple-DRM responses.    
    """
    def __init__(self):
        super().__init__()
        self._ebounds = None
        self._drms = []
    
    @property
    def detector(self):
        """(str): The name of the detector this response is associated with"""
        if len(self._drms) > 0:
            return self._drms[0].detector
        
    @property
    def ebounds(self):
        """(:class:`~.data_primitives.Ebounds`): The energy channel bounds"""
        return self._ebounds
    
    @property
    def num_chans(self):
        """(int): Number of energy channels"""
        if len(self._drms) > 0:
            return self._drms[0].num_chans

    @property
    def num_drms(self):
        """(int): Number of DRMs in the file"""
        return len(self._drms)

    @property
    def num_ebins(self):
        """(int): Number of photon energy bins"""
        if len(self._drms) > 0:
            return self._drms[0].num_ebins

    @property
    def tcent(self):
        """(np.array): The center times of the intervals over which the DRMs 
        are valid"""
        return (self.tstart + self.tstop) / 2.0

    @property
    def trigtime(self):
        """(float): The trigger time, if applicable"""
        return self._drms[0].trigtime

    @property
    def tstart(self):
        """(np.array): The start times of the intervals over which the DRMs 
        are valid"""
        return np.array([drm.tstart for drm in self._drms])

    @property
    def tstop(self):
        """(np.array): The end times of the intervals over which the DRMs 
        are valid"""
        return np.array([drm.tstop for drm in self._drms])

    def drm_index(self, time_range):
        """Return the indices of the DRMs that cover the requested time range.
        If the time range is before (after) the times covered by the DRMs in 
        the RSP2 object, then the first (last) index will be returned.
        
        Args:
            time_range (float, float): The time range
        
        Returns:        
            (np.array)
        """
        start, stop = self._assert_range(time_range)
        mask = (self.tstop > start) & (self.tstart < stop)
        
        if mask.sum() == 0:
            if stop < self.tstart[0]:
                return np.array([0], dtype=int)
            else:
                return np.array([self.num_drms-1], dtype=int)
        
        idx = np.arange(self.num_drms, dtype=int)
        return idx[mask]

    def extract_drm(self, index=None, atime=None):
        """Extract a single DRM from the Rsp2 object and return a single-DRM 
        :class:`Rsp` object. Either ``index`` or ``atime`` should be defined. 
        
        Args:
            index (int, optional): The DRM index to retrieve.
            atime (float, optional): The time corresponding to a DRM.  
                                     The nearest DRM to the time will be chosen.
        
        Returns:
            (:class:`Rsp`)
        """
        if (atime is None) and (index is None):
            raise RuntimeError('Either index or atime should be defined')
        if atime is not None:
            index = self.drm_index((atime, atime))[0]
        return self._drms[index]
    
    def interpolate(self, atime, **kwargs):
        """Interpolate over a multi-DRM object for a given time and return
        a new single-DRM RSP object. This function uses `scipy.interpolate.interp1d
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html>`_.
        
        If the given time is before (after) the times covered by the DRMs
        within the RSP object, then the first (last) DRM will be returned.
        
        Args:
            atime (float): The time corresponding to a DRM.  The nearest DRM to 
                           the time will be chosen.
            **kwargs: Keyword arguments to be passed to scipy.interpolate.interp1d
        
        Returns:
            (:class:`Rsp`)
        """
        if self.num_drms == 1:
            raise ValueError('Single DRM response.  Cannot interpolate')

        # create the interpolate and interpolate
        matrices = [rsp.drm.matrix for rsp in self._drms]
        matrix_interp = interp1d(self.tcent, np.array(matrices), axis=0,
                                 fill_value=(matrices[0], matrices[-1]),
                                 bounds_error=False, **kwargs)
        new_matrix = matrix_interp(atime)
        
        
        # create the object
        nearest_drm = self.nearest_drm(atime)
        new_drm = ResponseMatrix(new_matrix, 
                                 nearest_drm.drm.photon_bins.low_edges(),
                                 nearest_drm.drm.photon_bins.high_edges(),
                                 self.ebounds.low_edges(), 
                                 self.ebounds.high_edges())
        
        cls = type(nearest_drm)
        tstart = tstop = atime
        if nearest_drm.trigtime is not None:
            tstart += nearest_drm.trigtime
            tstop += nearest_drm.trigtime
        
        if nearest_drm.headers is None:
            headers = None
        else:
            headers = nearest_drm.headers.copy()
            
        rsp = cls.from_data(new_drm, start_time=tstart, stop_time=tstop,
                            trigger_time=nearest_drm.trigtime,
                            headers=headers)
        rsp._fchan = nearest_drm._fchan
        rsp._nchan = nearest_drm._nchan
        rsp._ngrp = nearest_drm._ngrp
        rsp._detector = self.detector
        return rsp

    def nearest_drm(self, atime):
        """Return a single-DRM Rsp object containing the nearest DRM to the 
        requested time.
        
        Args:
            atime (float): The requested time
        
        Returns: 
            (:class:`Rsp`)
        """
        idx = self.drm_index((atime, atime))[0]
        return self._drms[idx]

    def weighted(self, time_bins, interpolate=False, **kwargs):
        """Return a counts-weighted DRM.
        For long spectral accumulations, it is often appropriate to weight the
        responses by the count history. This is done by either selecting the
        DRM nearest to each time bin or interpolating the DRMs at each time 
        bin, weighting each DRM by the number of counts, summing the weighted
        DRMs, and normalizing by the total number of counts. A single-DRM Rsp
        object is returned.
        
        Args:
            time_bins (:class:`~.data_primitives.TimeBins`): 
                Time history count rate (lightcurve) data.
            interpolate (bool, optional): Set to True to interpolate at each
                                          time bin instead of selecting the 
                                          nearest DRM to each time bin.
                                          Default is False.
            **kwargs: Keyword arguments to be passed to scipy.interpolate.interp1d.
                      These will be ignored if ``interpolate=False``.
        
        Returns:        
            (:class:`Rsp`)
        """
        # only one DRM in the file; return it
        if self.num_drms == 1:
            return self[0]
        
        if not isinstance(time_bins, TimeBins):
            raise TypeError('time_bins must be of type TimeBins')
        
        # extract the corresponding DRM for each time bin and weight by
        # number of counts in the time bin
        bin_cents = time_bins.centroids
        counts = time_bins.counts
        matrix = np.zeros(self[0].drm.matrix.shape)
        for i in range(time_bins.size):
            if not interpolate:
                matrix += self.nearest_drm(bin_cents[i]).drm.matrix * counts[i]
            else:
                matrix += self.interpolate(bin_cents[i], **kwargs).drm.matrix * \
                          counts[i]
        matrix /= counts.sum()
        
        # create new DRM
        new_drm = ResponseMatrix(matrix, 
                                 self[0].drm.photon_bins.low_edges(),
                                 self[0].drm.photon_bins.high_edges(),
                                 self.ebounds.low_edges(), 
                                 self.ebounds.high_edges())
        
        # create new Response object
        cls = type(self[0])
        tstart, tstop = time_bins.range
        if self[0].trigtime is not None:
            tstart += self[0].trigtime
            tstop += self[0].trigtime
        
        if self[0].headers is None:
            headers = None
        else:
            headers = self[0].headers.copy()
            
        rsp = cls.from_data(new_drm, start_time=tstart, stop_time=tstop,
                            trigger_time=self[0].trigtime, headers=headers)
        
        return rsp
    
    def write(self, directory, filename=None, **kwargs):
        """Writes a multiple-DRM Rsp2 object to a FITS file.

        This method will not work unless the Rsp class for the responses 
        contained within the Rsp2 has been subclassed and the 
        :meth:`Rsp._build_hdulist` defined.
        
        Args:
            directory (str): The directory to write the file
            filename (str, optional): If set, will override the standardized name
            **kwargs (optional): keywords passed to 
                                 :meth:`astropy.io.fits.HDUList.write_to()`
        """
        # set filename
        if (self.filename is None) and (filename is None):
            raise NameError('Filename not set')
        
        if filename is None:
            filename = self.filename
        
        # update the creation time in the headers
        for rsp in self._drms:
            if rsp.headers is not None:
                rsp.headers.update()
                
        # build HDUs from each of the RSPs        
        rsp_hdus = []
        for rsp in self._drms:
            hdu = rsp._build_hdulist()
            if rsp.trigtime is not None:
                for h in hdu:
                    h.header['TSTART'] = rsp.tstart + rsp.trigtime
                    h.header['TSTOP'] = rsp.tstop + rsp.trigtime
            rsp_hdus.append(hdu)
                
        # use the primary HDU from the first RSP
        hdulist = fits.HDUList()
        hdulist.append(rsp_hdus[0][0])
        hdulist[0].header['FILENAME'] = filename
        hdulist[0].header['CREATOR'] = self._drms[0].headers.creator()[1]
        hdulist[0].header['DRM_NUM'] = self.num_drms
        
        # use the ebounds HDU from the first RSP
        hdulist.append(rsp_hdus[0][1])
        
        # all the specresp extensions
        hdulist.extend([hdu[2] for hdu in rsp_hdus])
        
        # write to file
        full_path = os.path.join(directory, filename)
        try:
            hdulist.writeto(full_path, checksum=True, **kwargs)
        except Exception as e: print(e)
    
    @classmethod
    def from_rsps(cls, rsp_list, filename=None):
        """Create a Rsp2 object from a list of Rsp objects.
        
        Args:
            rsp_list (list of :class:`Rsp`): The single-DRM responses
            filename (str, optional): The filename of the object
                 
        Returns:
            (:class:`Rsp2`)
        """
        obj = cls()
        obj._filename = filename
        
        # set data and ebounds
        try:
            if not all([isinstance(rsp, Rsp) for rsp in rsp_list]):
                raise TypeError('rsp_list must be a list of Rsp objects')
        except:
            raise TypeError('rsp_list must be a list of Rsp objects')
        
        obj._drms = rsp_list
        obj._ebounds = rsp_list[0].ebounds
                       
        return obj

    def _assert_range(self, valrange):
        assert valrange[0] <= valrange[1], \
            'Range must be in increasing order: (lo, hi)'
        return valrange

    def __getitem__(self, index):
        return self._drms[index]
    
    def __repr__(self):
        s = '<{0}: '.format(self.__class__.__name__)
        if self.filename is not None:
            s += '{};'.format(self.filename)
        if self.trigtime is not None:
            s += '\n trigger time: {};'.format(self.trigtime)
        s += ' {} DRMs;'.format(self.num_drms)
        s += '\n time range ({0}, {1})>'.format(self.tstart[0], self.tstop[-1])
        return s
    
