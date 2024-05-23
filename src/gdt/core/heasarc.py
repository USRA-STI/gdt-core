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

__all__ = ['ProgressMixin', 'Ftp', 'Http', 'BaseFinder', 'FtpFinder', 'FileDownloader', 'BrowseCatalog']


class ProgressMixin:

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


class BaseProtocol(AbstractContextManager, ABC, ProgressMixin):
    """A base class for the protocol used to access remote files.
    
    Note:
        This class should not be directly instantiated, but rather inherited.
        The inherited class should define methods called 
        ``_cd()``, ``_ls``, ``download``, and ``download_url``.
        
    Parameters:
        progress: The progress bar
    """

    def __init__(self, progress: Progress = None):
        self._progress = progress
        self._file_list = []

    @property
    def files(self):
        """(list of str): The list of files in the current directory"""
        return self._file_list

    @property
    def num_files(self):
        """(int): Number of files in the current directory"""
        return len(self._file_list)

    def cd(self, path: Union[str, Path]):
        """Change directory
        
        Args:
            path: The remote directory path
        """
        self._validate_connection()
        try:
            self._file_list = self.ls(path)
            self._cd(path)
        except:
            self._file_list = []
            raise ValueError('{} is not a valid path'.format(path))

    def ls(self, path: Union[str, Path]):
        """List the contents of a directory associated with
        a data set.
        
        Args:
            path: The remote directory path

        Returns:
            (list of str)
        """
        self._validate_connection()
        try:
            files = self._ls(path)
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

    def get(self, download_dir: Union[str, Path], files: List[str],
            verbose: bool = True) -> List[Path]:
        """Downloads a list of files from the current directory.
        This function also returns a list of the downloaded file paths.

        Args:
            download_dir (str, Path): The download directory location
            files (list of str): The list of files to download
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.        
        
        Returns:
            (list)
        """
        if not isinstance(files, list):
            raise ValueError("files argument must be a list.")

        filepaths = []
        for file in files:
            file_path = self.download(file, download_dir, verbose)
            filepaths.append(file_path)

        return filepaths

    def _reconnect(self):
        """Reconnect after connection loss. Default behavior is pass since
        not all inheriting classes need to reconnect."""
        pass

    def _validate_connection(self):
        """Validate that a connection remains open. Default behavior is pass
        since not all protocols maintain a continuous connection."""
        pass

    def __exit__(self, __exc_type: Optional[Type[BaseException]], __exc_value: Optional[BaseException],
                 __traceback: Optional[TracebackType]) -> Optional[bool]:
        pass

    @abstractmethod
    def download(self, file: str, dest_dir: Union[str, Path], verbose: bool = True):
        """Downloads a single file from the current directory.

        Args:
            file (str): The file name to download
            dest_dir (Path): The download file location
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.        
        
        Returns:
            (Path)
        """
        pass
 
    @abstractmethod
    def download_url(self, url: str, dest_dir: Union[str, Path], verbose: bool = True):
        """Download a file from a remote site

        Args:
            url (str): The url of a file to download
            dest_dir (Path): The directory where the file will be written
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (Path)
        """
        pass

    @abstractmethod
    def _cd(self, path: str):
        """Internal call to change directory

        Args:
            path: The remote directory path
        """
        pass

    @abstractmethod
    def _ls(self, path: str):
        """Internal call to list the contents of a directory associated with
        a data set.

        Args:
            path: The remote directory path

        Returns:
            (list of str)
        """
        pass


class Ftp(BaseProtocol):
    """A base class for FTP interactions with the HEASARC archive.
    
    Parameters:
        args: The set of parameters needed to define the data path
        host (str, optional): The host of the FTP archive
    """

    def __init__(self, host='heasarc.gsfc.nasa.gov', progress: Progress = None):
        super().__init__(progress)
        self._host = host
        self._ftp = None

        # If host is None, then let's not continue with the connection.
        if host is not None:
            self.connect(host=host)

    def _ls(self, path: Union[str, Path]):
        """List the directory contents of an FTP directory associated with
        a data set.
        
        Args:
            path: The remote directory path

        Returns:
            (list of str)
        """
        return self._ftp.nlst(path)

    def _cd(self, path):
        self._ftp.cwd(path)

    def _validate_connection(self):
        if self._ftp is None:
            raise ConnectionError('The connection is closed.')

    def download(self, file: str, dest_dir: Union[str, Path], verbose = True):
        """Downloads a single file from the current directory.

        Args:
            file (str): The file name to download
            file_path (Path): The download file location
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.        
        
        """
        # Make sure the dest_dir is a Path
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        file_path = dest_dir.joinpath(file)
 
        self._validate_connection()
        self._ftp.voidcmd('TYPE I')
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

        return file_path

    def download_url(self, url: str, dest_dir: Union[str, Path], verbose: bool = True):
        """Download a file from a remote site"""
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

        self.cd(str(url_path.parent))
        file_path = self.get(download_dir=dest_dir,
                             files=[url_path.name],
                             verbose=verbose)

        return file_path

    def connect(self, host: str = None):
        """Attempt a connection
        """
        if host is not None:
            self._host = host
        self._ftp = FTP_TLS(host=self._host)
        self._ftp.login()
        self._ftp.prot_p()

    def _reconnect(self):
        """Attempt a reconnect in case connection was lost
        """
        self._ftp.close()
        self.connect()

    def disconnect(self):
        """Politely disconnect from the FTP server."""
        if self._ftp is not None:
            self._ftp.quit()
            self._ftp.close()
        self._ftp = None
        self._file_list = []

    def pwd_r(self) -> str:
        """Retrieve the current directory."""
        self._validate_connection()
        if self._ftp is not None:
            return self._ftp.pwd()

    def __del__(self):
        if self._ftp is not None:
            self._ftp.close()

    def __exit__(self, __exc_type: Optional[Type[BaseException]], __exc_value: Optional[BaseException],
                 __traceback: Optional[TracebackType]) -> Optional[bool]:
        self.disconnect()
        return None

    def __repr__(self):
        return '<{0}: host {1}>'.format(self.__class__.__name__, self._host)


class Http(BaseProtocol):
    """A base class for HTTP/HTTPS interactions with the HEASARC archive.
    
    Parameters:
        args: The set of parameters needed to define the data path
        url (str, optional): The url of the HTTP/HTTPS archive
    """

    def __init__(self, url='https://heasarc.gsfc.nasa.gov/FTP',
                 start_key='<a href="', end_key='">', table_key='Parent Directory</a>',
                 progress: Progress = None, context: ssl.SSLContext = None):
        super().__init__(progress)
        self._url = url
        # keys are used to parse the HTTP/HTTPS file index
        self._start_key = start_key
        self._end_key = end_key
        self._table_key = table_key
        self._context = context
        # need to track current directory
        self._cwd = None

    def _cd(self, path: Union[str, Path]):
        """Change directory
        
        Args:
            path: The remote directory path
        """
        self._cwd = path

    def _ls(self, path: Union[str, Path]):
        """List the directory contents of an FTP directory associated with
        a data set.
        
        Args:
            path: The remote directory path

        Returns:
            (list of str)
        """
        files = []
        page = urlopen(self._url + path, context=self._context)
        table = page.read().decode("utf-8").split(self._table_key)[1]
        for line in table.split("\n"):
            if self._start_key in line:
                file = line.split(self._start_key)[1].split(self._end_key)[0]
                files.append(file)
        return files

    def download(self, file: str, dest_dir: Union[str, Path], verbose: bool = True):
        """Downloads a single file from the current directory.

        Args:
            file (str): The file name to download
            file_path (Path): The download file location
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.        
        
        """
        if self._cwd is None:
            raise ValueError("User must first cd() into a directory.")
        if self._url is None:
            raise ValueError("User must define base url at init.")

        file_path = self.download_url(
            self._url + self._cwd + file, dest_dir, verbose)

        return file_path

    def download_url(self, url: str, dest_dir: Union[str, Path], verbose: bool = True):
        """Download a file from a HTTP(S) site"""
        url_p = urlparse(url)

        # verify the URL is for HTTP(S)
        if url_p.scheme not in ('http', 'https'):
            raise ValueError('URL does not begin with http:// or https://')

        file = Path(url_p.path).name

        # Make sure the dest_dir is a Path
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        file_path = dest_dir.joinpath(file)

        response = urlopen(url, context=self._context)
        with file_path.open('wb') as fp:
            if verbose:
                total_size = int(response.headers["Content-Length"])
                if self._progress is None:
                    progress = self._create_progress()
                    progress.start()
                else:
                    progress = self._progress

                task_id = progress.add_task("download", filename=file, total=total_size)
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

        return file_path

    def __repr__(self):
        return '<{0}: url {1}>'.format(self.__class__.__name__, self._url)


class BaseFinder(AbstractContextManager, ABC):
    """A base class for the interface to the HEASARC FTP archive.

    Note:
        This class should not be directly instantiated, but rather inherited.
        The inherited class should define a method called
        ``_construct_path()`` that accepts ``*args``, which are the user-defined
        parameters required to define the data path, and the method should
        return the data path as a string

    Parameters:
        args: The set of parameters needed to define the data path
        protocol (str, optional): The connection protocol. Default is HTTPS.
        **kwargs: keyword arguments passed to protocol constructor
    """
    def __init__(self, *args, protocol='HTTPS', **kwargs):
       self.protocol = protocol
       if protocol in ['HTTP', 'HTTPS']:
           self._protocol = Http(**kwargs)
       elif protocol == 'FTP':
           self._protocol = Ftp(**kwargs)
       else:
           raise ValueError("Unrecognized connection protocol " + protocol)
       self._args = None

       self.cd(*args)

    @property
    def files(self):
        """(list of str): The list of files in the current directory"""
        return self._protocol._file_list

    @property
    def num_files(self):
        """(int): Number of files in the current directory"""
        return self._protocol._num_files

    def cd(self, *args):
        path = self._construct_path(*args)
        self._protocol.cd(path)

    def ls(self, *args):
        path = self._construct_path(*args)
        return self._protocol.ls(path)

    def filter(self, filetype, extension):
        """Filters the directory for the requested filetype and extension

        Args:
            filetype (str): The type of file, e.g. 'cspec'
            extension (str): The file extension, e.g. '.pha'

        Returns:
            (list)
        """
        return self._file_filter(self.files, filetype, extension)

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

    def __exit__(self, __exc_type: Optional[Type[BaseException]], __exc_value: Optional[BaseException],
                 __traceback: Optional[TracebackType]) -> Optional[bool]:
        pass


class FtpFinder(BaseFinder):
    def __init__(self, *args, host='heasarc.gsfc.nasa.gov', progress: Progress = None):
        super().__init__(*args, method='FTP', host=host, progress=progress)


class FileDownloader(ProgressMixin):
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
