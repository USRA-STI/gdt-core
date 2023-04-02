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
import astropy.io.fits as fits
from astropy.time import Time

from gdt.core import __version__
import copy

__all__ = ['Header', 'FileHeaders']


class Header(fits.Header):
    """A FITS header, subclassed from astropy.io.fits.Header.  
    This class should be further sub-classed with the following the class
    variables:
    
        *  name - The name of the extension
        *  keywords - A list of tuples that defined the keywords, their default
                      values, and associated comments
    
    The keyword list is treated as exhaustive on initialization.  In other words.
    the keyword list defines the only possible keywords that are allowed in the
    header.  The default value (not None), defines the data type of the keyword, 
    and is strictly typed (e.g. if the default value is an int, but a string 
    value not convertible to int is given, an exception will be raised).
    
    Once initialized, the class behaves as any other astropy.io.fits.Header, 
    however, new keywords cannot be added after initialization.
    """
    def __init__(self, *args, **kwargs):
        
        # an extension name must be set
        if not hasattr(self, 'name'):
            raise AttributeError("Header must have class attribute 'name' " \
                                 "defined")
        
        if not hasattr(self, 'keywords'):
            self.keywords = []
            
        self._kw_types = {kw[0]: type(kw[1]) for kw in self.keywords}
        
        # fill the header with the keywords
        super().__init__()
        for keyword in self.keywords:
            self.append(keyword)
        self.keywords = None
        
        for key, val in kwargs.items():
            self[key] = val
 
        if 'EXTNAME' in self.keys():
            self['EXTNAME'] = self.name

    def __setitem__(self, key, val):
        
        # pass-through for COMMENT
        if isinstance(key, tuple):
            super().__setitem__(key, val)
            return
        
        _keys = [k.lower() for k in self.keys()]
        if key.lower() not in _keys:
            if key.replace('_', '-').lower() not in _keys:
                raise KeyError('{} keyword does not exist'.format(key))
            else:
                key = key.replace('_', '-')
        the_type = self._kw_types[key.upper()]
        if val is not None:       
            try:
                val_type = the_type(val)
            except:
                raise TypeError('Value for {} is of incorrect type'.format(key)) 
        else:
            val_type = None
        super().__setitem__(key, val_type)
        
    @staticmethod
    def creator():
        """The Creator card"""
        return ('CREATOR','Gamma-ray Data Tools {}'.format(__version__), 
                'Software and version creating file') 
    

class FileHeaders():
    """A collection of FITS headers. This class should be further sub-classed 
    with the class variables:
    
        *  _header_templates - A list of :class:`Header`, where each header is
                               a template of values.
    
    Once initialized, each header can be accessed either by the extension name,
    e.g., 'PRIMARY' or by extension index.
    """
    def __init__(self):
        
        if not hasattr(self, '_header_templates'):
            raise AttributeError('FileHeaders must have class attribute '\
                                 '_header_templates')
        
        self._headers = {h.name: h for h in copy.deepcopy(self._header_templates)}
        self._header_templates = None       
        self.update()
    
    @property
    def num_headers(self):
        """(int): The number of headers"""
        return len(self._headers)
    
    def copy(self):
        """Return a copy of the FileHeaders.
        
        Returns:        
            (:class:`FileHeaders`)
        """
        return type(self).from_headers([self[i] for i in \
                                        range(self.num_headers)])
    
    def keys(self):
        """A list of the extension names
        
        Returns:        
            (list)
        """
        return list(self._headers.keys())

    def update(self):
        """Update the 'DATE' keywords to the current time.
        """
        date = Time.now().utc.isot
        for hdr in self._headers.values():
            try:
                hdr['DATE'] = date
            except:
                pass
                
    @classmethod
    def from_headers(cls, headers):
        """Create a FileHeaders object from a list of Headers.  The header types
        and keyword types must match those in ``_header_templates`` or an 
        exception is raised.
        
        Args:
            headers (list): A list of :class:`Headers`
                 
        Returns:
            (:class:`FileHeaders`)
        """
        obj = cls()
        num_headers = len(headers)
        if num_headers != len(obj._headers):
            raise ValueError('Incorrect number of headers for ' \
                             '{}'.format(cls.__name__))
        
        for i in range(num_headers):
            cidx = 0
            hidx = 0
            for key in obj[i].keys():
                if (key == 'COMMENT'):
                    obj[i][key][cidx] = headers[i][key][cidx]
                    cidx += 1
                elif (key == 'HISTORY'):
                    obj[i][key][hidx] = headers[i][key][hidx]
                    hidx += 1
                else:
                    obj[i][key] = headers[i][key]
        
        return obj
        
    @staticmethod
    def creator():
        """The Creator card"""
        return Header.creator()

    def __getitem__(self, key):
        if isinstance(key, int):
            num_keys = len(self.keys())
            if key < 0 or key > num_keys-1:
                raise IndexError('Out of range for {} headers'.format(num_keys))
            return self._headers[self.keys()[key]]
        else:
            if key not in list(self._headers.keys()):
                raise KeyError('{} header does not exist'.format(key))
            return self._headers[key]

    def __repr__(self):
        s = '<{0}: {1} headers>'.format(self.__class__.__name__, 
                                        self.num_headers)
        return s
