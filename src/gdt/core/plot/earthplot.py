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
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    
from .plot import GdtPlot, SAA, EarthLine, EarthPoints
from .defaults import *
from gdt.core.geomagnetic import *
from gdt.core.time import time_range

__all__ = ['EarthPlot']

class EarthPlot(GdtPlot):
    """Class for plotting the Earth, primarily related to a spacecraft's 
    position in orbit.
    
    Note:
        This class requires installation of Cartopy.
    
    Parameters:
        lat_range ((float, float), optional): 
            The latitude range of the plot. Default value is (-90, 90)
        lon_range ((float, float), optional):
            The longitude range of the plot. Default value is (-180, 180).
        saa (:class:`~gdt.core.geomagnetic.SouthAtlanticAnomaly`, optional): 
            If set, displays the SAA polygon.
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`        
    """
    def __init__(self, lat_range=None, lon_range=None, saa=None, **kwargs):
        if lat_range is None:
            lat_range = (-90.0, 90.0)
        if lon_range is None:
            lon_range = (-180.0, 180.0)

        # scale figure size by lat and lon ranges
        xsize = (lon_range[1] - lon_range[0]) / 30.0
        ysize = (lat_range[1] - lat_range[0]) / 30.0
        figsize = (xsize, ysize)

        super().__init__(figsize=figsize, **kwargs)
        self._saa = None
        self._orbit = None
        self._sc = None

        self._ax.get_xaxis().set_visible(False)
        self._ax.get_yaxis().set_visible(False)
        self._ax.set_visible(False)
        
        # initialize the map, mercator projection
        proj = ccrs.PlateCarree()
        self._m = self.fig.add_axes([0.05, 0.1, 0.9, 0.80], projection=proj)
        self._m.set_xlim(lon_range)
        self._m.set_ylim(lat_range)
        self._m.coastlines()
        
        gl = self._m.gridlines(crs=proj, draw_labels=True, color='black')
        gl.xlocator = mticker.FixedLocator(np.arange(-180., 181., 30.))
        gl.ylocator = mticker.FixedLocator(np.arange(-90., 91.0, 30.))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': PLOTFONTSIZE}
        gl.ylabel_style = {'size': PLOTFONTSIZE}
        gl.top_labels = False
        gl.right_labels = False

        if saa is not None:
            if not isinstance(saa, SouthAtlanticAnomaly):
                raise TypeError('saa must be of type SouthAtlanticAnomaly')
            self._saa = SAA(saa, proj, self._m, color='darkred', alpha=0.4)

    @property
    def geoaxes(self):
        """(cartopy.mpl.geoaxes.GeoAxesSubplot): The Cartopy axes"""
        return self._m

    @property
    def orbit(self):
        """(:class:`~gdt.plot.plot.EarthLine`): The orbital path"""
        return self._orbit

    @property
    def saa(self):
        """(:class:`~gdt.plot.plot.SAA`): The SAA polygon"""
        return self._saa

    @property
    def spacecraft(self):
        """(:class:`~gdt.plot.plot.EarthPoints`): A plot point representing the
        spacecraft position"""
        return self._sc

    def add_spacecraft_frame(self, frame, tstart=None, tstop=None, step=1.0, 
                             trigtime=None, icon=None, **kwargs):
        """Add a SpacecraftFrame object to plot the orbital path and 
        optionally the spacecraft position at a particular time.

        Args:
            frame (:class:`~gdt.core.coords.SpacecraftFrame`):
                The spacecraft frame containing position information
            tstart: (astropy.time.Time, optional): 
                The start time to plot the orbit.  If not set, then uses the 
                first ``frame.obstime``.
            tstop: (astropy.time.Time, astropy.units.Quantity, or float, optional):
                Either a time value signifying the end of the range, the 
                duration of the time range as a Quantity with appropriate unit, 
                or an int/float representing the number of seconds.  If not set,
                then uses the last ``frame.obstime``.
            step (datetime.TimeDelta, astropy.units.Quantity, or float, optional):
                Either a TimeDelta value, Quantity with appropriate unit, or an
                int/float representing the number of seconds. Default is 1.0
            trigtime (astropy.time.Time):
                Set to a particular time of interest to plot the location.
            icon (:class:`~gdt.plot.plot.EarthPoints`): 
                A subclass of EarthPoints that defines a specialized spacecraft
                icon.  If not set, will default to standard matplotlib markers.
            **kwargs: Arguments to be passed to 
                :class:`~gdt.plot.plot.EarthPoints` if ``icon`` is not set.
        """
        if tstart is None:
            tstart = frame.obstime[0]
        if tstop is None:
            tstop = frame.obstime[-1]
                
        # get latitude and longitude over the time range and produce the orbit
        times = time_range(tstart, tstop, step)
        frame_interp = frame.at(times)
        lat = frame_interp.earth_location.lat.value
        lon = frame_interp.earth_location.lon.value
        self._orbit = EarthLine(lat, lon, self._m, color='black',
                                alpha=0.4, lw=5.0)

        # if trigtime is set, produce spacecraft location
        if trigtime is not None:
            frame_interp = frame.at(trigtime)
            lat = frame_interp.earth_location.lat.value
            lon = frame_interp.earth_location.lon.value
            if icon is None:
                self._sc = EarthPoints(lat, lon, self._m, **kwargs)
            else:
                self._sc = icon(lat, lon, self._m, **kwargs)
            

    def standard_title(self):
        """Add a standard plot title containing orbital position of the 
        spacecraft.
        """
        if self.spacecraft is not None:
            coord = self.spacecraft.coordinates
            title = 'Latitude, East Longitude: ({0}, {1})\n'.format(*coord)
            #mcl = '{:.2f}'.format(float(self._trig_mcilwain))
            #title += 'McIlwain L: {}'.format(mcl)
            self._m.set_title(title)
