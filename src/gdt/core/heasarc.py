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
import socket
import ssl
import time
import shutil
from abc import ABC, abstractmethod
from contextlib import AbstractContextManager
from ftplib import FTP_TLS
from pathlib import Path
from types import TracebackType
from typing import List, Union, Type, Optional
from urllib.request import urlopen
from urllib.parse import urlparse
import numpy as np
import astropy.io.fits as fits

from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

__all__ = ['FtpFinder', 'BrowseCatalog', 'Http', 'Ftp', 'FileDownloader']


class BaseFinder(AbstractContextManager, ABC):

    def __init__(self, progress: Progress = None):
        self._progress = progress

    @staticmethod
    def _create_progress() -> Progress:
        """Creates a default progress object."""
        return Progress(
            TextColumn("[bold blue]{task.fields[filename]}", justify="right"),
            BarColumn(bar_width=None),
            "[progress.percentage]{task.percentage:>3.1f}%", "•",
            DownloadColumn(), "•",
            TransferSpeedColumn(), "•",
            TimeRemainingColumn(),
        )


class FtpFinder(BaseFinder):
    """A base class for the interface to the HEASARC FTP archive.
    
    Note:
        This class should not be directly instantiated, but rather inherited.
        The inherited class should define a method called 
        ``_construct_path()`` that accepts ``*args``, which are the user-defined
        parameters required to define the data path, and the method should
        return the data path as a string
        
    Parameters:
        args: The set of parameters needed to define the data path
        host (str, optional): The host of the FTP archive
    """

    def __init__(self, *args, host='heasarc.gsfc.nasa.gov', progress: Progress = None):
        super().__init__(progress)
        self._host = host
        self._args = None
        self._ftp = None
        self._file_list = []

        # If host is None, then let's not continue with the connection.
        if host is not None:
            self.connect(*args, host=host)
        else:
            if len(args) > 0:
                raise ValueError('*args were given while host was None')

    def __del__(self):
        if self._ftp is not None:
            self._ftp.close()

    def __validate_connection(self):
        if self._ftp is None:
            raise ConnectionError('The connection is closed.')

    @property
    def files(self):
        """(list of str): The list of files in the current directory"""
        return self._file_list

    @property
    def num_files(self):
        """(int): Number of files in the current directory"""
        return len(self._file_list)

    def cd(self, *args):
        """Change directory
        
        Args:
            args: The set of parameters needed to define the data path
        """
        self._args = self._validate(*args)

    def filter(self, filetype, extension):
        """Filters the directory for the requested filetype and extension
        
        Args:
            filetype (str): The type of file, e.g. 'cspec'
            extension (str): The file extension, e.g. '.pha'

        Returns:
            (list)
        """
        return self._file_filter(self.files, filetype, extension)

    def ls(self, *args):
        """List the directory contents of an FTP directory associated with
        a data set.
        
        Args:
            args: The set of parameters needed to define the data path

        Returns:
            (list of str)
        """
        self.__validate_connection()
        path = self._construct_path(*args)
        try:
            files = self._ftp.nlst(path)
        except AttributeError:
            print('Connection appears to have failed.  Attempting to reconnect...')
            try:
                self._reconnect()
                print('Reconnected.')
                return self.ls(id)
            except:
                raise RuntimeError('Failed to reconnect.')
        except:
            raise FileNotFoundError('{} does not exist'.format(path))
        return sorted([os.path.basename(f) for f in files])

    def connect(self, *args, host: str = None):
        """Attempt a connection
        """
        if host is not None:
            self._host = host
        self._ftp = FTP_TLS(host=self._host)
        self._ftp.login()
        self._ftp.prot_p()
        if len(args) > 0:
            self._args = self._validate(*args)

    @abstractmethod
    def _construct_path(self, *args) -> str:
        """This method needs to be defined by the inheriting class.  The method
        shall accept all user-defined parameters that are required to define 
        the data path and shall return the data path as a string.

        Args:
            args: The set of parameters needed to define the data path   
        
        Returns:
            (str)
        """
        pass

    def _file_filter(self, file_list, filetype, extension):
        """Filters the directory for the requested filetype and extension
        
        Args:
            file_list (list): A list of files
            filetype (str): The type of file, e.g. 'cspec'
            extension (str): The file extension, e.g. '.pha'

        Returns:
            list: The filtered file list
        """
        files = [f for f in file_list if
                 (filetype in f) & (f.endswith(extension))]

        return files

    def get(self, download_dir: Union[str, Path], files: List[str],
            verbose: bool = True) -> List[Path]:
        """Downloads a list of files from the current FTP directory.
        This function also returns a list of the downloaded file paths.

        Args:
            download_dir (str, Path): The download directory location
            files (list of str): The list of files to download
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.        
        
        Returns:
            (list)
        """
        self.__validate_connection()
        # convert download_dir to a Path and create the directory
        download_dir = Path(download_dir)
        download_dir.mkdir(parents=True, exist_ok=True)

        # download each file
        filepaths = []
        for file in files:
            # have to save in self because this can't be passed as an argument
            # in the callback

            # download file
            self._ftp.voidcmd('TYPE I')
            file_path = download_dir.joinpath(file)
            with file_path.open('wb') as fp:
                if verbose:
                    if self._progress is None:
                        progress = self._create_progress()
                        progress.start()
                    else:
                        progress = self._progress

                    task_id = progress.add_task("download", filename=file,
                                                total=self._ftp.size(file))

                    # the callback function
                    def write_func(data: bytes):
                        fp.write(data)
                        progress.update(task_id, advance=len(data))

                    self._ftp.retrbinary('RETR ' + file, callback=write_func)

                    # If progress was created locally then stop it here.
                    if self._progress is None:
                        progress.stop()
                else:
                    self._ftp.retrbinary('RETR ' + file, callback=fp.write)

            filepaths.append(file_path)

        return filepaths

    def _reconnect(self):
        """Attempt a reconnect in case connection was lost
        """
        self._ftp.close()
        self.connect()

    def _validate(self, *args):
        """Validate arguments by constructing the FTP path and attempting to
        change to that directory"""
        self.__validate_connection()
        try:
            self._file_list = self.ls(*args)
            self._ftp.cwd(self._construct_path(*args))
            return args
        except:
            self._file_list = []
            raise ValueError('{} are not valid arguments'.format(args))

    def disconnect(self):
        """Politely disconnect from the FTP server."""
        if self._ftp is not None:
            self._ftp.quit()
            self._ftp.close()
        self._ftp = None
        self._file_list = []

    def pwd_r(self) -> str:
        """Retrieve the current directory."""
        self.__validate_connection()
        if self._ftp is not None:
            return self._ftp.pwd()

    def __exit__(self, __exc_type: Optional[Type[BaseException]], __exc_value: Optional[BaseException],
                 __traceback: Optional[TracebackType]) -> Optional[bool]:
        self.disconnect()
        return None

    def __repr__(self):
        args = ', '.join([str(arg) for arg in self._args]) \
            if self._args is not None else ''

        return '<{0}: {1}>'.format(self.__class__.__name__, args)


class Ftp(FtpFinder):

    def _construct_path(self, *args) -> str:
        return str(args[0])

    def download_url(self, url: str, dest_dir: Union[str, Path], verbose: bool = True):
        """Download a file from a FTP site"""
        url_p = urlparse(url)

        # verify the URL is for FTP
        if url_p.scheme != 'ftp':
            raise ValueError('URL does not begin with ftp://')

        url_path = Path(url_p.path)

        if self._host is None:
            self.connect(host=url_p.hostname)
        elif self._host != url_p.hostname:
            self.disconnect()
            self.connect(host=url_p.hostname)

        self.cd(url_path.parent)
        self.get(download_dir=dest_dir, files=[url_path.name], verbose=verbose)


class Http(BaseFinder):

    def download_url(self, url: str, dest_dir: Union[str, Path], verbose: bool = True):
        """Download a file from a HTTP(S) site"""
        url_p = urlparse(url)

        # verify the URL is for HTTP(S)
        if url_p.scheme not in ('http', 'https'):
            raise ValueError('URL does not begin with http:// or https://')

        url_path = Path(url_p.path)

        # Make sure the dest_dir is a Path
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        file = dest_dir.joinpath(url_path.name)

        response = urlopen(url)
        with file.open('wb') as fp:
            if verbose:
                total_size = int(response.headers["Content-Length"])
                if self._progress is None:
                    progress = self._create_progress()
                    progress.start()
                else:
                    progress = self._progress

                task_id = progress.add_task("download", filename=file.name, total=total_size)
                while True:
                    data = response.read(32768)
                    if not data:
                        break
                    fp.write(data)
                    progress.update(task_id, advance=len(data))

                # If this is a locally created progress, then let's stop it.
                if self._progress is None:
                    progress.stop()
            else:
                shutil.copyfileobj(response, fp)

    def __exit__(self, __exc_type: Optional[Type[BaseException]], __exc_value: Optional[BaseException],
                 __traceback: Optional[TracebackType]) -> Optional[bool]:
        pass

class FileDownloader(BaseFinder):
    """Used to download a list of files given as a URL."""

    def __init__(self):
        super().__init__(progress=self._create_progress())
        self._progress.start()
        self._http = Http(progress=self._progress)
        self._ftp = Ftp(host=None, progress=self._progress)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._ftp.disconnect()
        self._progress.stop()

    def download_url(self, url: str, dest_dir: Union[str, Path], verbose: bool = True):
        """Download a file from a URL."""

        if url.startswith('ftp'):
            self._ftp.download_url(url, dest_dir, verbose)
        elif url.startswith('http'):
            self._http.download_url(url, dest_dir, verbose)
        else:
            raise ValueError('url must begin with ftp://, http://, or https://')

    def bulk_download(self, urls: List[str], dest_dir: Union[str, Path], verbose: bool = True):
        """Download files from a list of URLs."""
        for url in urls:
            self.download_url(url, dest_dir, verbose)



class BrowseCatalog:
    """A class that interfaces with the HEASARC Browse API.  This can be
    called directly, but is primarily intended as a base class.
    
    This class makes a query to HEASARC's w3query.pl perl script in 
    BATCHRETRIEVALCATALOG mode.  All fields and rows are retrieved so that
    this class, on instantiation, contains the full set of catalog data. 
    Any queries based on row or columns selections/slices are then done locally,
    instead of making repeated requests to the HEASARC. A cached copy of the 
    catalog is saved locally so future instantiations can utilize the local copy
    instead of querying HEASARC.
    
    Parameters:
        cache_path (str): The path where the cached catalog will live.
        cached (bool, optional): Set to True to read from the cached file
                                 instead of querying HEASARC. Default is False.
        table (str, optional): The name of the table to be passed to the 
                               w3query.pl script.
        verbose (bool, optional): Default is True    
    """

    def __init__(self, cache_path, table=None, verbose=True, cached=False):

        if not os.path.exists(cache_path):
            os.makedirs(cache_path)
        self._cache_path = cache_path

        self._table = str(table)
        self._cached = bool(cached)
        self._verbose = bool(verbose)

        # just load the cached file
        if self._cached:
            self._header, self._data = self._read_cache()
            return

        # build the URL query
        host = 'https://heasarc.gsfc.nasa.gov'
        script = 'db-perl/W3Browse/w3query.pl'
        query = 'tablehead=name=BATCHRETRIEVALCATALOG_2.0+'
        # Retrieve all fields, all rows, and return in FITS format 
        query += '{}&Fields=All&displaymode=FitsDisplay&ResultMax=0'.format(table)

        if self._table is not None:
            self._is_connected(host)
            self._download_table(host + '/' + script + '?' + query)
            self._header, self._data = self._read_cache()

    @property
    def columns(self):
        """(np.array): The names of the columns available in the table"""
        return self._data.dtype.names

    @property
    def num_cols(self):
        """(int): The total number of columns (fields) in the data table"""
        return len(self.columns)

    @property
    def num_rows(self):
        """(int): The total number of rows in the data table"""
        return self._data.size

    def column_range(self, column):
        """Return the data range for a given column, in the form of (low, high).
        
        Args:
            column (str): The column name

        Returns:
            (tuple)
        """
        col_copy = np.copy(self._data[column])
        col_copy.sort()
        return col_copy[0], col_copy[-1]

    def get_table(self, columns=None):
        """Return the table data as a numpy record array.
        
        Args:
            columns (list of str, optional): The columns to return. If omitted, 
                                             returns all columns.

        Returns:
            (np.recarray)
        """
        if columns is None:
            columns = self.columns

        arrays = [self._data[col] for col in columns]
        return np.rec.fromarrays(arrays, names=','.join(columns).upper())

    def slice(self, column, lo=None, hi=None):
        """Perform row slices of the data table based on a conditional of a
        single column. Returns a new BrowseCatalog object.
        
        Args:
            column (str): The column name
            lo (optional): The minimum (inclusive) value of the slice. If not 
                           set, uses the lowest range of the data in the column.
            hi (optional): The maximum (inclusive) value of the slice. If not 
                           set, uses the highest range of the data in the column.

        Returns:
            (:class:`BrowseCatalog`)
        """
        col = self._data[column]
        if lo is None:
            lo, _ = self.column_range(column)
        if hi is None:
            _, hi = self.column_range(column)
        mask = (col >= lo) & (col <= hi)

        # create a new object and fill it with the sliced data
        obj = type(self)(self._cache_path, cached=True)
        obj._table = self._table
        obj._header = np.copy(self._header)
        obj._data = self._data[mask]
        return obj

    def slices(self, columns):
        """Perform row slices of the data table based on a conditional of 
        multiple columns. Returns a new BrowseCatalog object.
        
        Args:
            columns (list of tuples):
                A list of tuples, where each tuple is (column, lo, hi).  The 
                'column' is the column name, 'lo' is the lowest bounding value, 
                and 'hi' is the highest bouding value.  If no low or high 
                bounding is desired, set to None. See :meth:`slice()` for more 
                info.

        Returns:
            (:class:`BrowseCatalog`)
        """
        numcols = len(columns)
        obj = self
        for i in range(numcols):
            obj = obj.slice(columns[i][0], lo=columns[i][1], hi=columns[i][2])
        return obj

    def _download_table(self, url):
        """Downloads the table from HEASARC
        """
        # secure connection
        if self._verbose:
            print('Sending request and awaiting response from HEASARC...')
        t0 = time.time()

        context = ssl._create_unverified_context()
        page = urlopen(url, context=context)

        if self._verbose:
            print(f'Downloading {self._table} from HEASARC via w3query.pl...')

        with open(os.path.join(self._cache_path, self._table + '.fit'),
                  'wb') as f:
            f.write(page.read())

        if self._verbose:
            print('Finished in {} s'.format(int(time.time() - t0)))

    def _is_connected(self, host):
        """Test the connection to the host to determine if it is reachable"""
        try:
            sock = socket.create_connection((host.split('/')[-1], 80))
            sock.close()
        except OSError:
            raise OSError("Either you are not connected to the internet or "
                          "{0} is down.".format(host))
        return True

    def _read_cache(self):
        """Read the cached catalog FITS file"""
        out_file = os.path.join(self._cache_path, self._table + '.fit')
        with fits.open(out_file) as hdulist:
            return hdulist[1].header, hdulist[1].data

    def __repr__(self):
        return '<{0}: {1} columns, {2} rows>'.format(self.__class__.__name__,
                                                     self.num_cols,
                                                     self.num_rows)
