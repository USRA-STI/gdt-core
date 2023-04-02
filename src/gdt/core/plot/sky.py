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
from astropy.coordinates import get_sun, SkyCoord
from astropy.time import Time
from gdt.core.detector import Detectors
from .plot import DetectorPointing, GalacticPlane, SkyHeatmap, SkyPolygon
from .plot import PlotElementCollection as Collection
from .plot import GdtPlot, SkyLine, SkyCircle, Sun
from .lib import *

__all__ = ['SkyPlot', 'EquatorialPlot', 'GalacticPlot', 'SpacecraftPlot']


def get_lonlat(skycoord):
    """Given a SkyCoord, retrieve the longitude/latitude values in that order.
    Because different frames can have different names for the lon/lat variables,
    there is not a consistent frame-agnostic way of retrieving these values.
    
    Args:
        skycoord (astropy.coordinates.SkyCoord): The skycoord
    
    Returns:
        (astropy.coordinates.angle.Longitude, 
         astropy.coordinates.angle.Latitude)  
    """
    names = skycoord.get_representation_component_names()
    lon = None
    lat = None
    for k, v in names.items():
        if v == 'lon':
            lon = getattr(skycoord, k)
        elif v == 'lat':
            lat = getattr(skycoord, k)
        else:
            pass
    return (lon, lat)


class SkyPlot(GdtPlot):
    """Base class for an all-sky plot.

    Parameters:
        projection (str, optional): The projection of the map. Current 
                                    compatible projections: 'aitoff', 'hammer',
                                    'lambert', 'mollweide', and 'polar'. 
                                    Default is 'mollweide'.
        flipped (bool, optional): 
            If True, the longitudinal axis is flipped, following 
            astronomical convention. Default is True.
        xticks_res (float, optional): The resolution, in degrees, of the 
                                      longitudinal tick marks. Default is 30.
        yticks_res (float, optional): The resolution, in degrees, of the 
                                      latitudinal tick marks. Default is 15.
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`
    """
    _background = 'antiquewhite'
    _textcolor = 'black'
    _canvascolor = 'white'
    _fontsize = 10
    _frame = None
    _astropy_frame = None

    def __init__(self, projection='mollweide', flipped=True, xticks_res=30, 
        yticks_res=15, **kwargs):
        super().__init__(figsize=(10, 5), projection=projection,
                         **kwargs)

        # set up the plot background color and the sky grid
        self._figure.set_facecolor(self._canvascolor)
        self._ax.set_facecolor(self._background)
        self._ax.grid(True, linewidth=0.5)

        # create the axes tick labels
        self._longitude_axis(flipped, xticks_res)
        self._latitude_axis(yticks_res)

        self._flipped = flipped
        self._sun = None
        self._earth = None
        self._detectors = Collection()
        self._galactic_plane = None
        self._posterior = None
        self._eff_area = None
        self._clevels = Collection()

    @property
    def detectors(self):
        """(:class:`~gdt.plot.plot.PlotElementCollection` of \
        :class:`~gdt.plot.plot.DetectorPointing`): The collection of detector 
        plot elements"""
        return self._detectors

    @property
    def earth(self):
        """(:class:`~gdt.plot.plot.SkyCircle`): The Earth plot element"""
        return self._earth
    
    @property
    def effective_area(self):
        """(:class:`~gdt.plot.plot.SkyHeatmap`): The effective area plot element"""
        return self._eff_area
           
    @property
    def fontsize(self):
        """(int): The font size of the text labels."""
        return self._fontsize
    @fontsize.setter
    def fontsize(self, size):
        self._ax.set_yticklabels(self._ytick_labels, fontsize=size,
                                 color=self._textcolor)
        self._ax.set_xticklabels(self._xtick_labels, fontsize=size,
                                 color=self._textcolor)
        self._fontsize = size

    @property
    def galactic_plane(self):
        """(:class:`~gdt.plot.plot.GalacticPlane`): The galactic plane plot 
        element"""
        return self._galactic_plane

    @property
    def loc_contours(self):
        """(:class:`~gdt.plot.plot.PlotElementCollection` of \
        :class:`~gdt.plot.plot.SkyLine` or :class:`~gdt.plot.plot.SkyPolygon`):
        The localization contour plot elements"""
        return self._clevels

    @property
    def loc_posterior(self):
        """(:class:`~gdt.plot.plot.SkyHeatmap`): The localization gradient plot 
        element"""
        return self._posterior

    @property
    def sun(self):
        """(:class:`~gdt.plot.plot.Sun`): The Sun plot element"""
        return self._sun

    @property
    def text_color(self):
        """(str): The color of the text labels"""
        return self._textcolor
    @text_color.setter
    def text_color(self, color):
        self._ax.set_yticklabels(self._ytick_labels, fontsize=self._fontsize,
                                 color=color)
        self._ax.set_xticklabels(self._xtick_labels, fontsize=self._fontsize,
                                 color=color)
        self._textcolor = color

    def add_frame(self, frame, trigtime=None, detectors='all', earth=True,
                    sun=True, galactic_plane=True):
        """Add a SpacecraftFrame object to plot the location of the Earth, Sun, 
        and detector pointings.

        Args:
            frame (:class:`gdt.core.coords.SpacecraftFrame`): The spacecraft
                frame containing position and orientation information
            trigtime (astropy.time.Time, optional): 
                The time of interest. This must be set if ``frame`` contains 
                more than one frame with a corresponding time.  If ``frame`` is
                of size=1, then the corresponding time for that frame will be
                used.                               
            detectors (list or "all"): A list of detectors or "all" to plot the 
                                       pointings on the sky.
            earth (bool, optional): If True, plot the Earth. Default is True.
            sun (bool, optional): If True, plot the Sun. Default is True.
            galactic_plane (bool, optional):
                If True, plot the Galactic plane. Default is True.
        """        
        # if SpacecraftFrame is of only one size, use that trigtime
        if frame.ndim == 0:
            trigtime = frame.obstime
        else:
            if trigtime is None:
                raise ValueError('trigtime must be set or a SpacecraftFrame ' \
                                 'corresponding to a single time must be used')
        
        if trigtime is not None:
            if not isinstance(trigtime, Time):
                raise TypeError('trigtime must either be an astropy Time object')
            
            if frame.ndim != 0:
                frame_interp = frame.at(trigtime)
            else:
                frame_interp = frame
            
            if sun:
                self.plot_sun(trigtime)
            if earth:
                self.plot_earth(frame_interp.geocenter, 
                                frame_interp.earth_angular_radius)

            if detectors == 'all':
                detectors = [det.name for det in frame.detectors]
            else:
                if isinstance(detectors, str):
                    detectors = [detectors]
            det_objs = [frame.detectors.from_str(detector) \
                        for detector in detectors]
            for det in det_objs:
                det_coord = det.skycoord(frame).transform_to('icrs')
                self.plot_detector(det_coord, det.name)

        if galactic_plane:
            self.plot_galactic_plane()

    def add_effective_area(self, hpx, frame=None, sun=False, earth=False, 
                           detectors=[], galactic_plane=False):
        """Add a HealPixEffectiveArea object to plot the effective area. 
        Optionally add a SpacecraftFrame object to plot the location of the 
        Earth, Sun, galactic plane, and detector pointings.

        Args:
            hpx (:class:`~gdt.core.healpix.HealPixEffectiveArea`): 
                The HEALPix object
            frame (:class:`gdt.core.coords.SpacecraftFrame`, optional): 
                The spacecraft frame containing position and orientation 
                information
            detectors ('all' or list):
                A list of detectors or "all" to plot the pointings on the sky
            earth (bool, optional): If True, plot the Earth. Default is True.
            sun (bool, optional): If True, plot the Sun. Default is True.
            galactic_plane (bool, optional):
                If True, plot the Galactic plane. Default is True.
        """
        # determine what the resolution of the sky grid should be based on the
        # resolution of the healpix
        approx_res = np.sqrt(hpx.pixel_area)
        numpts_az = int(np.floor(0.5*360.0/approx_res))
        numpts_zen = int(np.floor(0.5*180.0/approx_res))

        eff_arr, az_arr, zen_arr = hpx.make_grid(numpts_az=numpts_az,
                                                 numpts_zen=numpts_zen)
        
        if frame is not None:
            x, y = np.meshgrid(az_arr, 90.0-zen_arr)
            coords = SkyCoord(x.flatten(), y.flatten(), frame=frame, unit='deg').gcrs
            ra = coords.ra.value.reshape(x.shape)
            dec = coords.dec.value.reshape(y.shape)
            self._eff_area = self.plot_heatmap(eff_arr, ra, dec)
        else:
            self._eff_area = self.plot_heatmap(eff_arr, az_arr, 90.0-zen_arr)
    
        if frame is not None:
            if sun:
                self.plot_sun(frame.obstime)
            
            if earth:
                self.plot_earth(frame.geocenter, frame.earth_angular_radius)
            
            if detectors == 'all':
                detectors = [det.name for det in frame.detectors]
            else:
                if isinstance(detectors, str):
                    detectors = [detectors]
            det_objs = [frame.detectors.from_str(detector) \
                        for detector in detectors]
            for det in det_objs:
                det_coord = det.skycoord(frame)
                self.plot_detector(det_coord, det.name)
            
            if galactic_plane:
                self.plot_galactic_plane()

    def add_localization(self, hpx, gradient=True, clevels=None, sun=True,
                         earth=True, detectors='all', galactic_plane=True):
        """Add a HealPixLocalization object to plot a localization and 
        optionally plot the location of the Earth, Sun, and detector pointings.

        Args:
            hpx (:class:`~gdt.core.healpix.HealPixLocalization`): 
                The HEALPix object
            gradient (bool, optional): 
                If True, plots the posterior as a color gradient. If False, 
                plot the posterior as color-filled confidence regions.
            clevels (list of float, optional):
                The confidence levels to plot contours. Default plots at the 
                1, 2, and 3 sigma level.
            detectors ('all' or list):
                A list of detectors or "all" to plot the pointings on the sky
            earth (bool, optional): If True, plot the Earth. Default is True.
            sun (bool, optional): If True, plot the Sun. Default is True.
            galactic_plane (bool, optional):
                If True, plot the Galactic plane. Default is True.
        
        Note:
            Setting `gradient=False` when plotting an annulus may produce 
            unexpected results at this time.  It is suggested to use 
            `gradient=True` for plotting annuli maps.
        """
        if clevels is None:
            clevels = [0.997, 0.955, 0.687]
        
        if detectors == 'all':
            detectors = [det.name for det in hpx.frame.detectors]
        else:
            if isinstance(detectors, str):
                detectors = [detectors]
        det_objs = [hpx.frame.detectors.from_str(detector) \
                    for detector in detectors]
        
        # determine what the resolution of the sky grid should be based on the
        # resolution of the healpix
        approx_res = np.sqrt(hpx.pixel_area)
        numpts_ra = int(np.floor(0.5*360.0/approx_res))
        numpts_dec = int(np.floor(0.5*180.0/approx_res))
        
        # plot the gradient
        if gradient:
            prob_arr, ra_arr, dec_arr = hpx.prob_array(numpts_ra=numpts_ra,
                                                       numpts_dec=numpts_dec)
            self._posterior = self.plot_heatmap(prob_arr, ra_arr, dec_arr)
        
        # plot the confidence levels.
        for clevel in clevels:
            paths = hpx.confidence_region_path(clevel, numpts_ra=numpts_ra, 
                                               numpts_dec=numpts_dec)
            
            lons = []
            lats = []
            for path in paths:
                coord = SkyCoord(path[:,0], path[:,1], frame='gcrs', unit='deg')
                lon, lat = get_lonlat(coord.transform_to(self._astropy_frame))
                lons.append(lon)
                lats.append(lat)
            
            # if we plotted the gradient, then only plot the unfilled contours,
            # otherwise, plot the stacked, filled contours
            numpaths = len(paths)
            if gradient:
                for i in range(numpaths):
                    contour = SkyLine(lons[i].value, lats[i].value, self.ax,
                                      frame=self._frame, color='black',
                                      alpha=0.7, linewidth=2,
                                      flipped=self._flipped)
                    self._clevels.include(contour, str(clevel) + '_' + str(i))
            else:
                for i in range(numpaths):
                    contour = SkyPolygon(lons[i].value, lats[i].value, self.ax,
                                         frame=self._frame, color='purple',
                                         alpha=None, face_alpha=0.3, 
                                         flipped=self._flipped)
                    self._clevels.include(contour, str(clevel) + '_' + str(i))
        
        # important for polar plot...maybe other projections, too
        try:
            self._ax.set_rlim(-np.pi/2.0, np.pi/2.0)
        except:
            pass
        
        # plot sun
        if sun:
            try:
                self.plot_sun(hpx.frame.obstime)
            except:
                pass
        
        # plot earth
        if earth:
            try:
                self.plot_earth(hpx.frame.geocenter, 
                                hpx.frame.earth_angular_radius)
            except:
                pass
        
        # plot detector pointings
        for det in det_objs:
            det_coord = det.skycoord(hpx.frame)
            self.plot_detector(det_coord, det.name)
        
        # plot galactic plane
        if galactic_plane:
            self.plot_galactic_plane()
        
    def plot_detector(self, det_coord, det, radius=10.0, **kwargs):
        """Plot a detector pointing.

        Args:
            det_coord (astropy.coordinates.SkyCoord): The detector coordinate
            det (str): The detector name
            radius (float, optional): The radius of pointing, in degrees. 
                                      Default is 10.0
            **kwargs: Options to pass to :class:`~gdt.plot.plot.SkyCircle`
        """
        x, y = get_lonlat(det_coord.transform_to(self._astropy_frame))
        pointing = DetectorPointing(x.value, y.value, radius, det, 
                                    self.ax, frame=self._frame, 
                                    flipped=self._flipped, **kwargs)
        self._detectors.include(pointing, det)

    def plot_earth(self, geo_coord, radius, **kwargs):
        """Plot the Earth.

        Args:
            geo_coord (astropy.coordinates.SkyCoord): The coordinates of the 
                                                      geocenter.
            radius (astropy.quantity.Quantity): The angular radius of the Earth
            **kwargs: Options to pass to :class:`~gdt.plot.plot.SkyCircle`
        """
        geo_coord = SkyCoord(geo_coord.ra, geo_coord.dec, frame='gcrs')
        lon, lat = get_lonlat(geo_coord.transform_to(self._astropy_frame))
        radius = radius.to('deg').value
        self._earth = SkyCircle(lon.value, lat.value, radius, self.ax, 
                                flipped=self._flipped, frame=self._frame, 
                                color='deepskyblue', alpha=None, 
                                face_alpha=0.25, edge_alpha=0.50, **kwargs)

    def plot_galactic_plane(self):
        """Plot the Galactic plane.
        """
        self._galactic_plane = GalacticPlane(self.ax, flipped=self._flipped,
                                             frame=self._frame)

    def plot_heatmap(self, heatmap, lon_array, lat_array, **kwargs):
        """Plot a heatmap on the sky.

        Args:
            heatmap (np.array): A 2D array of values
            ra_array (np.array): The array of longitudinal gridpoints
            dec_array (np.array): The array of latitudinal gridpoints
            **kwargs: Options to pass to :class:`~gdt.plot.plot.SkyHeatmap`
        """
        if lon_array.ndim == 1 and lat_array.ndim == 1:
            x, y = np.meshgrid(lon_array, lat_array)
        else:
            x = lon_array
            y = lat_array
        coords = SkyCoord(x.flatten(), y.flatten(), frame='gcrs', unit='deg')
        
        try:
            coords = coords.transform_to(self._astropy_frame)
        except:
            pass
        
        lon, lat = get_lonlat(coords)
        lon = lon.reshape(x.shape)
        lat = lat.reshape(y.shape)
        heatmap = SkyHeatmap(lon.value, lat.value, heatmap, self.ax, 
                             frame=self._frame, flipped=self._flipped, **kwargs)
        return heatmap

    def plot_sun(self, time, **kwargs):
        """Plot the sun.

        Args:
            time (astropy.time.Time): The time at which to plot the sun
            **kwargs: Options to pass to :class:`~gdt.plot.plot.Sun`
        """
        sun_coord = get_sun(time)
        sun_coord = SkyCoord(sun_coord.ra, sun_coord.dec, 
                             frame='gcrs').transform_to(self._astropy_frame)
        lon, lat = get_lonlat(sun_coord)
        self._sun = Sun(lon.value, lat.value, self.ax, 
                        flipped=self._flipped, frame=self._frame, **kwargs)

    def _longitude_axis(self, flipped, xticks_res):
        # longitude labels
        # these have to be shifted on the plot because matplotlib natively
        # goes from -180 to +180
        ticks = np.linspace(-np.pi, np.pi, int(360.0/xticks_res)+1)
        self._ax.set_xticks(ticks)

        # important for polar plot...maybe other projections, too
        try:
            self._ax.set_thetalim(-np.pi, np.pi)
        except:
            pass
        
        tick_labels = np.linspace(0, 360, int(360.0/xticks_res)+1)
        tick_labels += self._x_start
        tick_labels = np.remainder(tick_labels, 360)
        if flipped:
            tick_labels = tick_labels[::-1]
        # format the tick labels with degrees
        self._xtick_labels = [str(int(t)) + '$^\circ$' for t in tick_labels]
        # remove label to avoid overlap with latitude labels
        self._xtick_labels[0] = ''
        self._ax.set_xticklabels(self._xtick_labels, fontsize=self._fontsize,
                                 color=self._textcolor)

    def _latitude_axis(self, yticks_res):
        # latitude labels
        # matplotlib natively plots from -90 to +90 from bottom to top. 
        # this is fine for equatorial coordinates, but we have to shift if 
        # we are plotting in spacecraft coordinates
        ticks = np.linspace(-np.pi/2.0, np.pi/2.0, int(180.0/yticks_res)+1)
        self._ax.set_yticks(ticks)

        tick_labels = np.linspace(-90, 90, int(180.0/yticks_res)+1)
        tick_labels += self._y_start        
        if np.sign(self._y_start) == -1:
            tick_labels *= -1
        self._ytick_labels = [str(int(t)) + '$^\circ$' for t in tick_labels]
        self._ax.set_yticklabels(self._ytick_labels, fontsize=self._fontsize,
                                 color=self._textcolor)


class EquatorialPlot(SkyPlot):
    """Plotting the sky in Equatorial (GCRS) coordinates.

    Parameters:
        projection (str, optional): The projection of the map. Current 
                                    compatible projections: 'aitoff', 'hammer',
                                    'lambert', 'mollweide', and 'polar'. 
                                    Default is 'mollweide'.
        flipped (bool, optional): 
            If True, the longitudinal axis is flipped, following 
            astronomical convention. Default is True.
        xticks_res (float, optional): The resolution, in degrees, of the 
                                      longitudinal tick marks. Default is 30.
        yticks_res (float, optional): The resolution, in degrees, of the 
                                      latitudinal tick marks. Default is 15.
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`
    """
    _x_start = 0
    _y_start = 0
    _frame = 'equatorial'
    _astropy_frame = 'gcrs'

    def add_effective_area(self, hpx, frame=None, **kwargs):
        if frame is None:
            raise ValueError('frame must be set to plot in Equatorial')
        super().add_effective_area(hpx, frame=frame, **kwargs)
 

class GalacticPlot(SkyPlot):
    """Plotting the sky in Galactic coordinates.

    Parameters:
        projection (str, optional): The projection of the map. Current 
                                    compatible projections: 'aitoff', 'hammer',
                                    'lambert', 'mollweide', and 'polar'. 
                                    Default is 'mollweide'.
        flipped (bool, optional): 
            If True, the longitudinal axis is flipped, following 
            astronomical convention. Default is True.
        xticks_res (float, optional): The resolution, in degrees, of the 
                                      longitudinal tick marks. Default is 30.
        yticks_res (float, optional): The resolution, in degrees, of the 
                                      latitudinal tick marks. Default is 15.
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`
    """
    _x_start = 180
    _y_start = 0
    _frame = 'galactic'
    _astropy_frame = 'galactic'

    def add_effective_area(self, hpx, frame=None, **kwargs):
        if frame is None:
            raise ValueError('frame must be set to plot in Galactic')
        super().add_effective_area(hpx, frame=frame, **kwargs)


class SpacecraftPlot(SkyPlot):
    """Class for plotting the sky in Spacecraft coordinates.

    Parameters:
        zenith (bool, optional): Set to True to plot the latitudinal coordinates
                                 relative to spacecraft zenith, otherwise plots
                                 relative to spacecraft elevation. Default is True
        projection (str, optional): The projection of the map. Current 
                                    compatible projections: 'aitoff', 'hammer',
                                    'lambert', 'mollweide', and 'polar'. 
                                    Default is 'mollweide'.
        flipped (bool, optional): 
            If True, the longitudinal axis is flipped, following 
            astronomical convention. Default is True.
        xticks_res (float, optional): The resolution, in degrees, of the 
                                      longitudinal tick marks. Default is 30.
        yticks_res (float, optional): The resolution, in degrees, of the 
                                      latitudinal tick marks. Default is 15.
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`
    """
    _x_start = 180
    _y_start = -90
    _frame = 'spacecraft'
    
    def __init__(self, zenith=True, **kwargs):
        if zenith:
            self._y_start = -90
        else:
            self._y_start = 0
        super().__init__(flipped=False, **kwargs)
        
    def add_frame(self, frame, **kwargs):
        self._astropy_frame = frame
        super().add_frame(frame, **kwargs)

    def add_effective_area(self, hpx, frame=None, **kwargs):
        self._astropy_frame = frame
        super().add_effective_area(hpx, frame=frame, **kwargs)
    
    def add_localization(self, hpx, **kwargs):
        self._astropy_frame = hpx.frame
        super().add_localization(hpx, **kwargs)
    
    def plot_galactic_plane(self):
        """Plot the Galactic plane.
        """
        self._galactic_plane = GalacticPlane(self.ax, flipped=self._flipped,
                                             frame=self._astropy_frame)

