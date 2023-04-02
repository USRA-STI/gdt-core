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
import astropy.io.fits as fits

from .data_primitives import Ebounds, Gti, EnergyBins
from .file import FitsFileContextManager
from .headers import Header, FileHeaders
from gdt.core.background.primitives import BackgroundSpectrum

__all__ = ['Pha', 'Bak']

# header information

_chantype_card = ('CHANTYPE', 'PHA', 'No corrections have been applied')
_detchans_card = ('DETCHANS', 0, 'Total number of channels in each rate')
_extname_card = ('EXTNAME', '', 'name of this binary table extension')
_filter_card = ('FILTER', '', 'The instrument filter in use (if any)')
_hduclass_card = ('HDUCLASS', 'OGIP', 'Conforms to OGIP standard indicated in HDUCLAS1')
_hduvers_card = ('HDUVERS', '1.2.1', 'Version of HDUCLAS1 format in use')
_trigtime_card = ('TRIGTIME', 0.0, 'Trigger time')

class PrimaryHeader(Header):
    name = 'PRIMARY'
    keywords = [Header.creator(), 
                ('DATE', '', 'file creation date (YYYY-MM-DDThh:mm:ss UT)'),
                ('TSTART', 0.0, 'Observation start time'),
                ('TSTOP', 0.0, 'Observation stop time'),
                _trigtime_card,
                ('FILENAME', '', 'Name of this file')]

class EboundsHeader(Header):
    name = 'EBOUNDS'
    keywords = [_extname_card, _hduclass_card,
                ('HDUCLAS1', 'RESPONSE', 'These are typically found in RMF ' \
                                         'files'),
                ('HDUCLAS2', 'EBOUNDS', 'From CAL/GEN/92-002'), 
                _hduvers_card, _chantype_card, _filter_card, _detchans_card,
                ('CH2E_VER', '', 'Channel to energy conversion scheme used'),
                ('GAIN_COR', 0.0, 'Gain correction factor applied to energy ' \
                                  'edges')]

class PhaSpectrumHeader(Header):
    name = 'SPECTRUM'
    keywords = [_extname_card, _filter_card,
            ('AREASCAL', 1., 'No special scaling of effective area by channel'),
            ('BACKFILE', '', 'Name of corresponding background file (if any)'),
            ('BACKSCAL', 1., 'background file scaling factor'),
            ('CORRFILE', '', 'Name of corresponding correction file (if any)'),
            ('CORRSCAL', 1., 'Correction scaling file'),
            ('RESPFILE', '', 'Name of corresponding RMF file (if any)'),
            ('ANCRFILE', '', 'Name of corresponding ARF file (if any)'),
            ('SYS_ERR', 0., 'Systematic errors'),
            ('POISERR', True, 'Assume Poisson Errors'),
            ('GROUPING', 0, 'No special grouping has been applied'),
            _hduclass_card,
            ('HDUCLAS1', 'SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)'),
            ('HDUCLAS2', 'TOTAL', 'Indicates gross data (source + background)'),
            ('HDUCLAS3', 'COUNT', 'Indicates data stored as counts'),
            ('HDUCLAS4', 'TYPEI', 'Indicates PHA Type I file format'),
            _hduvers_card, _chantype_card, _detchans_card, 
            ('EXPOSURE', 0.0, 'Accumulation time - deadtime')]

class BakSpectrumHeader(Header):
    name = 'SPECTRUM'
    keywords = [_extname_card, _filter_card,
            ('AREASCAL', 1., 'No special scaling of effective area by channel'),
            ('BACKFILE', '', 'Name of corresponding background file (if any)'),
            ('BACKSCAL', 1., 'background file scaling factor'),
            ('CORRFILE', '', 'Name of corresponding correction file (if any)'),
            ('CORRSCAL', 1., 'Correction scaling file'),
            ('RESPFILE', '', 'Name of corresponding RMF file (if any)'),
            ('ANCRFILE', '', 'Name of corresponding ARF file (if any)'),
            ('SYS_ERR', 0., 'Systematic errors'),
            ('POISERR', True, 'Assume Poisson Errors'),
            ('GROUPING', 0, 'No special grouping has been applied'),
            _hduclass_card,
            ('HDUCLAS1', 'SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)'),
            ('HDUCLAS2', 'BKG', 'Background PHA Spectrum'),
            ('HDUCLAS3', 'RATE', 'PHA data stored as rates'),
            ('HDUCLAS4', 'TYPEI', 'Indicates PHA Type I file format'),
            _hduvers_card, _chantype_card, _detchans_card, 
            ('EXPOSURE', 0.0, 'Accumulation time - deadtime')]

class GtiHeader(Header):
    name = 'GTI'
    keywords = [_extname_card, _hduclass_card,
                ('HDUCLAS1', 'GTI', 'Indicates good time intervals'),
                _hduvers_card]

class PhaHeaders(FileHeaders):
    _header_templates = [PrimaryHeader(), EboundsHeader(), PhaSpectrumHeader(), 
                         GtiHeader()]

class BakHeaders(FileHeaders):
    _header_templates = [PrimaryHeader(), EboundsHeader(), BakSpectrumHeader(), 
                         GtiHeader()]


class Pha(FitsFileContextManager):
    """PHA class for count spectra.
    """
    def __init__(self):
        super().__init__()
        self._ebounds = None
        self._data = None
        self._gti = None
        self._trigtime = None
        self._channel_mask = None

    @property
    def channel_mask(self):
        """(np.array): A Boolean array representing the valid channels"""
        if self._channel_mask is None:
            return np.ones(self.num_chans, dtype=bool)
        else:
            return self._channel_mask

    @property
    def data(self):
        """(:class:`~.data_primitives.EnergyBins`): The PHA data"""
        return self._data

    @property
    def ebounds(self):
        """(:class:`~.data_primitives.Ebounds`): The energy bounds"""
        return self._ebounds

    @property
    def energy_range(self):
        """(float, float): The energy range of the spectrum"""
        if self._data is not None:
            return self._data.range

    @property
    def exposure(self):
        """(float): The exposure of the PHA data"""
        return self._data.exposure[0]

    @property
    def gti(self):
        """(:class:`~.data_primitives.Gti`): The good time intervals"""
        return self._gti

    @property
    def num_chans(self):
        """(int): The number of energy channels"""
        return self._data.size

    @property
    def tcent(self):
        """(float): The center time of the data"""
        return sum(self.time_range) / 2.0

    @property
    def time_range(self):
        """(float, float): The time range of the spectrum"""
        if self._gti is not None:
            return self._gti.range

    @property
    def trigtime(self):
        """(float): The trigger time of the data, if available"""
        return self._trigtime

    @property
    def valid_channels(self):
        """(np.array): The channel indices that are valid"""
        return np.arange(self.num_chans, dtype=int)[self.channel_mask]

    def rebin_energy(self, method, *args, energy_range=(None, None), **kwargs):
        """Rebin the PHA in energy given a rebinning method. 

        Args:
            method (<function>): The rebinning function
            *args: Arguments to be passed to the rebinning function
            energy_range ((float, float), optional): 
                The starting and ending energy to rebin.  If omitted, uses the 
                full range of data.  Setting start or end to ``None`` will use 
                the data from the beginning or end of the data, respectively.
        Returns        
            (:class:`Pha`)
        """
        emin, emax = self._assert_range(energy_range)
        data = self.data.rebin(method, *args, emin=emin, emax=emax)
        headers = self._build_headers(self.trigtime, *self.time_range, 
                                      data.size)
        
        pha = self.from_data(data, gti=self.gti, trigger_time=self.trigtime, 
                             headers=headers, **kwargs)
        return pha

    def slice_energy(self, energy_ranges, **kwargs):
        """Slice the PHA by one or more energy ranges. Produces a new 
        PHA object.

        Args:
            energy_ranges ([(float, float), ...]): 
                The energy ranges to slice the data to.
        
        Returns:        
            (:class:`Pha`)
        """
        energy_ranges = self._assert_range_list(energy_ranges)
        data = [self.data.slice(*self._assert_range(energy_range)) \
                for energy_range in energy_ranges]
        data = EnergyBins.merge(data)

        headers = self._build_headers(self.trigtime, *self.time_range, 
                                      data.size)
        
        pha = self.from_data(data, gti=self.gti, trigger_time=self.trigtime, 
                               headers=headers, **kwargs)
        
        return pha

    @classmethod
    def from_data(cls, data, gti=None, trigger_time=None, filename=None,
                  headers=None, channel_mask=None, header_type=PhaHeaders,
                  **kwargs):
        """Create a PHA object from an 
        :class:`~.data_primitives.EnergyBins` object.

        Args:
            data (:class:`~.data_primitives.EnergyBins`): 
                 The PHA count spectrum data
            gti (:class:`~.data_primitives.Gti`), optional): 
                The good time intervals of the pectrum data.  If omitted, then 
                assumes the range (0, exposure).
            trigger_time (float, optional): The trigger time, if applicable. 
                                            If provided, the data times will be 
                                            shifted relative to the trigger time.
                                            Default is zero.
            headers (:class:`~.headers.FileHeaders`): The file headers
            channel_mask (np.array(dtype=bool)): 
                A boolean array representing the valid channels. If omitted, 
                assumes all non-zero count channels are valid.
            header_type (:class:`~.headers.FileHeaders`): 
                Default file header class. Only used if ``headers`` is not 
                defined
        
        Returns:
            (:class:`Pha`)
        """
        obj = cls()
        obj._filename = filename
        
        # set data and ebounds
        if not isinstance(data, EnergyBins):
            raise TypeError('data must be of type EnergyBins')
        obj._data = data
        obj._ebounds = Ebounds.from_bounds(data.lo_edges, data.hi_edges)
        
        # set GTI
        if gti is not None:
            if not isinstance(gti, Gti):
                raise TypeError('gti must be of type Gti')
        else:
            gti = Gti.from_list([(0.0, data.exposure[0])])
        obj._gti = gti

        # update times to be relative to ...if trigtime is set
        if trigger_time is not None:
            if trigger_time < 0.0:
                raise ValueError('trigger_time must be non-negative')
            obj._trigtime = trigger_time

        # set headers
        if headers is not None:
            if not isinstance(headers, FileHeaders):
                raise TypeError('headers must be of type FileHeaders')
            obj._headers = headers
        else:
            obj._headers = header_type()   
            tstart, tstop = obj._gti.range
            obj._headers['PRIMARY']['TSTART'] = tstart
            obj._headers['PRIMARY']['TSTOP'] = tstop
            obj._headers['EBOUNDS']['DETCHANS'] = data.size  
            obj._headers['SPECTRUM']['DETCHANS'] = data.size
        obj._headers['PRIMARY']['FILENAME'] = filename
        obj._headers['PRIMARY']['TRIGTIME'] = trigger_time
        obj._headers['SPECTRUM']['EXPOSURE'] = obj._data.exposure[0]
        
        # set the channel mask
        # if no channel mask is given, assume zero-count channels are bad
        if channel_mask is None:
            channel_mask = np.zeros(data.size, dtype=bool)
            channel_mask[data.counts > 0] = True
        try:
            iter(channel_mask)
            channel_mask = np.asarray(channel_mask).flatten().astype(bool)
        except:
            raise TypeError('channel_mask must be a Boolean array')
        if channel_mask.size != obj._data.size:
            raise ValueError('channel_mask must be the same size as the ' \
                             'number of data bins')
        obj._channel_mask = channel_mask
        
        return obj

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a PHA FITS file and return the PHA object
        
        If this class is inherited, this method may be over-written if a 
        non-standard file is being parsed, or if there is extra header 
        information/data that needs to be stored.
                
        Args:
            file_path (str): The file path of the FITS file
        
        Returns:        
            (:class:`Pha`)
        """
        obj = super().open(file_path, **kwargs)
        trigtime = None

        # get the headers
        hdrs = [hdu.header for hdu in obj.hdulist]
        if hdrs[2]['HDUCLAS2'] == 'TOTAL':
            headers = PhaHeaders.from_headers(hdrs)
        elif hdrs[2]['HDUCLAS2'] == 'BKG':
            headers = BakHeaders.from_headers(hdrs)

        if 'TRIGTIME' in hdrs[0].keys():
            trigtime = float(headers['PRIMARY']['TRIGTIME'])
            
        # data
        exposure = headers['SPECTRUM']['EXPOSURE']            
        if headers['SPECTRUM']['HDUCLAS3'] == 'COUNT':
            data = EnergyBins(obj.column(2, 'COUNTS'), obj.column(1, 'E_MIN'),
                              obj.column(1, 'E_MAX'), exposure)
        elif headers['SPECTRUM']['HDUCLAS3'] == 'RATE':
            data = BackgroundSpectrum(obj.column(2, 'RATES'), 
                                      obj.column(2, 'STAT_ERR'),
                                      obj.column(1, 'E_MIN'), 
                                      obj.column(1, 'E_MAX'), exposure)
        else:
            raise ValueError('Unable to determin PHA type. Looking for' \
                             '"COUNT" or "RATE" value for HDUCLAS3')
                
        # Quality flag indicates which channels are valid.
        # Valid: Quality = 0
        # Invalid: Quality = 1 
        if 'QUALITY' in obj.get_column_names(2):
            channel_mask = (obj.column(2, 'QUALITY') == 0)
        else:
            channel_mask = None
            
        # GTI
        gti_start = obj.column(3, 'START')
        gti_stop = obj.column(3, 'STOP')
        if trigtime is not None:
            gti_start -= trigtime
            gti_stop -= trigtime
        gti = Gti.from_bounds(gti_start, gti_stop)
            
        obj._hdulist.close()
        
        # create object
        obj = cls.from_data(data, gti=gti, trigger_time=trigtime, 
                            filename=obj.filename, 
                            headers=headers, channel_mask=channel_mask)
        
        
        return obj

    def _assert_range(self, valrange):
        if valrange[0] is None and valrange[1] is None:
            return valrange
        assert valrange[0] <= valrange[1], \
            'Range must be in increasing order: (lo, hi)'
        return valrange

    def _assert_range_list(self, range_list):
        try:
            iter(range_list[0])
        except:
            range_list = [tuple(range_list)]
        range_list = [self._assert_range(r) for r in range_list]
        return range_list
    
    def _build_hdulist(self):
        """This builds the HDU list for the FITS file.  
        
        If this class is inherited, this method may be over-written if a 
        non-standard file is being written, or if there is extra header 
        information/data that needs to be written.
        
        This method should construct each HDU (PRIMARY, EBOUNDS, SPECTRUM, GTI, 
        etc.) containing the respective header and data. The HDUs should then 
        be inserted into a HDUList and that list returned
        
        Returns:
            (:class:`astropy.io.fits.HDUList`)
        """
        # create FITS and primary header
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU(header=self.headers['PRIMARY'])
        for key, val in self.headers['PRIMARY'].items():
            primary_hdu.header[key] = val
        hdulist.append(primary_hdu)
        
        # the ebounds extension
        ebounds_hdu = self._ebounds_table()
        hdulist.append(ebounds_hdu)
        
        # the spectrum extension
        spectrum_hdu = self._spectrum_table()
        hdulist.append(spectrum_hdu)        
        
        # the GTI extension
        gti_hdu = self._gti_table()
        hdulist.append(gti_hdu)
        
        return hdulist

    def _build_headers(self, trigtime, tstart, tstop, num_chans):
        """This builds the headers for the FITS file.  
        
        If this class is inherited, this method may be over-written if there is 
        extra header information that needs to be written.
        
        This method should construct the headers from the minimum required 
        arguments and additional keywords and return a 
        :class:`~.headers.FileHeaders` object.
        
        Args:
            trigtime (float or None): The trigger time.  Set to None if no
                                      trigger time.
            tstart (float): The start time
            tstop (float): The stop time
            num_chans (int): Number of detector energy channels
        
        Returns:
            (:class:`~.headers.FileHeaders`)
        """
        headers = self.headers.copy()
        headers['PRIMARY']['TRIGTIME'] = trigtime
        headers['PRIMARY']['TSTART'] = tstart
        headers['PRIMARY']['TSTOP'] = tstop
        headers['EBOUNDS']['DETCHANS'] = num_chans
        headers['SPECTRUM']['DETCHANS'] = num_chans        
        return headers

    def _ebounds_table(self):
        chan_col = fits.Column(name='CHANNEL', format='1I', 
                               array=np.arange(self.num_chans, dtype=int))
        emin_col = fits.Column(name='E_MIN', format='1E', unit='keV', 
                               array=self.ebounds.low_edges())
        emax_col = fits.Column(name='E_MAX', format='1E', unit='keV', 
                               array=self.ebounds.high_edges())
        
        hdu = fits.BinTableHDU.from_columns([chan_col, emin_col, emax_col], 
                                            header=self.headers['EBOUNDS'])
        for key, val in self.headers['EBOUNDS'].items():
            hdu.header[key] = val

        return hdu

    def _spectrum_table(self):
        chan_col = fits.Column(name='CHANNEL', format='1I', 
                               array=np.arange(self.num_chans, dtype=int))
        counts_col = fits.Column(name='COUNTS', format='J', bzero=32768, 
                                 bscale=1, unit='count', array=self.data.counts)
        qual_col = fits.Column(name='QUALITY', format='1I', 
                               array=(~self.channel_mask).astype(int))
        
        hdu = fits.BinTableHDU.from_columns([chan_col, counts_col, qual_col], 
                                            header=self.headers['SPECTRUM'])
        for key, val in self.headers['SPECTRUM'].items():
            hdu.header[key] = val
        hdu.header.comments['TZERO2'] = 'offset for unsigned integers'
        hdu.header.comments['TSCAL2'] = 'data are not scaled'
        return hdu

    def _gti_table(self):
        tstart = np.array(self.gti.low_edges())
        tstop = np.array(self.gti.high_edges())
        if self.trigtime is not None:
            tstart += self.trigtime
            tstop += self.trigtime

        start_col = fits.Column(name='START', format='1D', unit='s', 
                                bzero=self.trigtime, array=tstart)
        stop_col = fits.Column(name='STOP', format='1D', unit='s', 
                                bzero=self.trigtime, array=tstop)
        hdu = fits.BinTableHDU.from_columns([start_col, stop_col], 
                                            header=self.headers['GTI'])
        
        for key, val in self.headers['GTI'].items():
            hdu.header[key] = val        
        hdu.header.comments['TZERO1'] = 'Offset, equal to TRIGTIME'
        hdu.header.comments['TZERO2'] = 'Offset, equal to TRIGTIME'
        return hdu

    def __repr__(self):
        s = '<{0}: '.format(self.__class__.__name__)
        if self.filename is not None:
            s += '{};'.format(self.filename)
        if self.trigtime is not None:
            s += '\n trigger time: {};'.format(self.trigtime)
        s += '\n time range {};'.format(self.time_range)
        s += '\n energy range {}>'.format(self.energy_range)
        return s


class Bak(Pha):
    """Class for a PHA background spectrum.    
    """
    @classmethod
    def from_data(cls, data, gti=None, trigger_time=None, filename=None,
                  headers=None, channel_mask=None, **kwargs):
        """Create a BAK object from an 
        :class:`~..background.background.BackgroundSpectrum` object.

        Args:
            data (:class:`~..background.background.BackgroundSpectrum`): 
                 The background spectrum data
            gti (:class:`~.data_primitives.Gti`), optional): 
                The good time intervals of the pectrum data.  If omitted, then 
                assumes the range (0, exposure).
            trigger_time (float, optional): The trigger time, if applicable. 
                                            If provided, the data times will be 
                                            shifted relative to the trigger time. Default is zero.
            headers (:class:`~.headers.FileHeaders`): The file headers
            channel_mask (np.array(dtype=bool)): 
                A boolean array representing the valid channels. If omitted, 
                assumes all non-zero count channels are valid.
        
        Returns:
            (:class:`Bak`)
        """
        return super().from_data(data, gti=gti, trigger_time=trigger_time, 
                                 filename=filename, headers=headers, 
                                 channel_mask=channel_mask, 
                                 header_type=BakHeaders, **kwargs)

    def rebin_energy(self, method, *args, emin=None, emax=None):
        """Not Implemented"""
        raise NotImplementedError('Function not available for BAK objects')

    def slice_energy(self, method, *args, emin=None, emax=None):
        """Not Implemented"""
        raise NotImplementedError('Function not available for BAK objects')

    def _spectrum_table(self):
        chan_col = fits.Column(name='CHANNEL', format='1I', 
                               array=np.arange(self.num_chans, dtype=int))
        rates_col = fits.Column(name='RATES', format='1D', unit='count/s', 
                                array=self.data.rates)
        staterr_col = fits.Column(name='STAT_ERR', format='1D', unit='count/s', 
                                  array=self.data.rate_uncertainty)
        
        hdu = fits.BinTableHDU.from_columns([chan_col, rates_col, staterr_col], 
                                            header=self.headers['SPECTRUM'])
        for key, val in self.headers['SPECTRUM'].items():
            hdu.header[key] = val
        return hdu

