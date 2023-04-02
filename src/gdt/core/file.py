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
from contextlib import AbstractContextManager
from pathlib import Path
from typing import List, Union

import numpy as np
from astropy.io import fits

__all__ = ['FileContextManager', 'FitsFileContextManager']


class FileContextManager(AbstractContextManager):
    def __init__(self, file_path: Union[str, Path], mode: str = 'r'):
        path = Path(file_path)
        self.file_obj = path.open(mode)

    def close(self):
        if self.file_obj is not None:
            self.file_obj.close()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __repr__(self):
        return f'<{self.__class__.__name__}(filename="{self.file_obj.name}", mode="{self.file_obj.mode}") ' \
               f'at {hex(id(self))}>'


class FitsFileContextManager(AbstractContextManager):
    """A context manager for FITS files.  Includes some convenience functions.
    """

    def __init__(self):
        self._filename = None
        self._hdulist = None
        self._headers = None

    @property
    def filename(self):
        """(str): The filename"""
        return self._filename

    @property
    def hdulist(self):
        """(astropy.io.fits.hdu.HDUList): The list of Header Data Units"""
        if self._hdulist is not None:
            return self._hdulist
        else:
            return self._build_hdulist()

    @property
    def headers(self):
        """(:class:`~gdt.core.headers.FileHeaders`): The headers"""
        return self._headers

    @property
    def num_hdus(self):
        """(int): The number of HDUs"""
        return len(self.hdulist)

    def close(self):
        """Close the file"""
        if self._hdulist is not None:
            self._hdulist.close()

    def column(self, hdu_num: int, col_name: str) -> np.array:
        """Return a column from an HDU as an array.
    
        Args:
            hdu_num (int): The HDU number
            col_name (str): The name of the column
    
        Returns:
            (np.array)
        """
        return np.array(self.hdulist[hdu_num].data[col_name])

    def columns_as_array(self, hdu_num: int, col_names: List[str],
                         dtype: np.dtype = None) -> np.array:
        """Return a list of columns from an HDU as an array.
    
        Args:
            hdu_num (int): The HDU number
            col_names (list of str): The names of the columns
            dtype (np.dtype, optional): The custom dtype of the output array
    
        Returns:
            (np.array)
        """
        return np.array([self.column(hdu_num, x) for x in col_names],
                        dtype=dtype).T

    def get_column_names(self, hdu_num: int):
        """Get the column names in a HDU.  Returns empty if there is no data
        extension in the HDU.
    
        Args:
            hdu_num (int): The HDU number
    
        Returns:
            (tuple)
        """
        try:
            return self.hdulist[hdu_num].data.dtype.names
        except:
            return ()

    @classmethod
    def open(cls, file_path: Union[str, Path], mode: str = 'readonly', memmap: bool = None):
        """Open a FITS file.
            
        Args:
            file_path (str): The file path
            mode (str): The file access mode
            memmap (bool): If True, memory map when reading the file
        
        Returns:
            (:class:`FitsFileContextManager`)
        """
        path = Path(file_path)
        obj = cls()
        obj._hdulist = fits.open(path, mode=mode, memmap=memmap)
        obj._filename = path.name
        return obj

    def write(self, directory: Union[str, Path], filename: str = None, **kwargs):
        """Write the file to disk.
        
        Args:
            directory (str): The directory to write the file.
            filename (str, optional): The filename.  If omitted, attempts to use
                                      the :attr:`~FitsFileContextManager.filename`
                                      if set.
        """
        if (self.filename is None) and (filename is None):
            raise NameError('Filename not set')

        dir_path = Path(directory)

        if filename is None:
            filename = self.filename

        try:
            self.headers['PRIMARY']['FILENAME'] = filename
        except:
            pass

        try:
            self.headers['PRIMARY']['CREATOR'] = self.headers.creator()[1]
        except:
            pass

        # update the creation time in the headers
        if self.headers is not None:
            self.headers.update()

        # write to file
        full_path = dir_path / filename
        try:
            self.hdulist.writeto(full_path, checksum=True, **kwargs)
        except Exception as e:
            print(e)

    def _build_hdulist(self):
        """This builds the HDU list for the FITS file.  This method needs
        to be specified in the inherited class.  The method should construct
        each HDU (PRIMARY, EBOUNDS, SPECTRUM, GTI, etc.) containing the 
        respective header and data. The HDUs should then be inserted into a
        HDUList and that list returned
        
        Returns:
            (:class:`astropy.io.fits.HDUList`)
        """
        raise NotImplementedError

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __repr__(self):
        return f'<{self.__class__.__name__}(filename="{self.filename}") at {hex(id(self))}>'

    def _repr_html_(self):
        s = f'<p>&lt{self.__class__.__name__}(filename=<b>"{self.filename}"</b>) at {hex(id(self))}&gt</p>'
        s += '<table>'
        s += '<tr><th>No.</th><th>Name</th><th>Ver</th><th>Type</th><th>Cards</th><th>Dimensions</th></tr>'
        for row in self._hdulist.info(False):
            s += f'<tr><td>{row[0]}</td><td>{row[1]}</td><td>{row[2]}</td><td>{row[3]}</td><td>{row[4]}</td>' \
                 f'<td>{row[5]}</td></tr>'
        s += '</table>'
        return s
