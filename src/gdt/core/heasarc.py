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

#__all__ = ['FtpFinder', 'BrowseCatalog', 'Http', 'Ftp', 'FileDownloader']
__all__ = ['ProgressMixin', 'Ftp', 'Http']


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

    def __init__(self, *args, url='https://heasarc.gsfc.nasa.gov/FTP',
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

