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
import re
import copy
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import colorConverter, PowerNorm
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
from .lib import *
from gdt.core.collection import DataCollection

__all__ = ['DetectorPointing', 'EarthLine', 'EarthPoints', 'EffectiveArea',
           'GalacticPlane', 'GdtCmap', 'GdtPlot', 'Heatmap', 'Histo', 
           'HistoErrorbars', 'HistoFilled', 'LightcurveBackground',
           'ModelData', 'ModelSamples', 'PlotElement', 'PlotElementCollection',
           'SAA', 'SkyAnnulus', 'SkyHeatmap', 'SkyCircle', 'SkyLine', 
           'SkyPoints', 'SkyPolygon', 'SpectrumBackground', 'Sun']

class GdtPlot():
    """Base class for plots in the GDT.
    
    Note:
        This class is not intended to be instantiated directly by the user, 
        rather it is inherited by different plot classes.
    
    Parameters:
        figsize ((float, float), optional): The figure size (width, height)
        dpi (int, optional): The resolution of the plot. Default is 100
        canvas (Canvas Backend object, optional): 
            If interfacing with a backend (e.g. Tk), pass the relavant 
            Canvas Backend object and the master widget window. If not set, 
            uses the default the matplotlib backend.
        interactive (bool, optional): If True, then enables interactive plotting 
                                     (python environment is available while plot 
                                     is displayed). Default is False. This is 
                                     overriden if canvas is set.
        ax (:class:`matplotlib.axes`, optional): 
            An existing axes object to plot to.  If not set, will create a new 
            axes object.
        **kwargs: Additional options for `Figure.add_subplot 
            <https://matplotlib.org/3.2.0/api/_as_gen/matplotlib.figure.Figure.html#matplotlib.figure.Figure.add_subplot>`_.    
    """
    def __init__(self, figsize=None, dpi=100, canvas=None, interactive=False,
                 ax=None, **kwargs):
        self._canvas = None
        self._figure = None
        
        if ax is not None:
            self._ax = ax
        else:
            if canvas is not None:
                if figsize is None:
                    figsize = (5, 4)
                self._figure = Figure(figsize=figsize, dpi=dpi)
                # if tkagg, remember to assign to frame
                self._canvas = canvas[0](self._figure, canvas[1])
            else:
                # just use pyplot, but disable blocking to enable simultaneous
                # access to the python command line and dynamic updating
                if figsize is None:
                    figsize = (7.7, 4.7)
                self._figure = plt.figure(figsize=figsize, dpi=dpi)
                self._canvas = self._figure.canvas
                if interactive:
                    plt.ion()
            self._ax = self._figure.add_subplot(111, **kwargs)


        self._format_coordinates()

    @property
    def ax(self):
        """(:class:`matplotlib.axes`): The matplotlib axes object for the plot"""
        return self._ax

    @property
    def canvas(self):
        """(Canvas Backend object): The plotting canvas, if set upon 
        initialization."""
        return self._canvas

    @property
    def fig(self):
        """(:class:`matplotlib.figure`): The matplotlib figure object"""
        return self._figure

    @property
    def xlim(self):
        """(float, float): The plotting range of the x axis"""
        return self._ax.get_xlim()
    @xlim.setter
    def xlim(self, val):
        self._ax.set_xlim(val)

    @property
    def xscale(self):
        """(str): The scale of the x axis, either 'linear' or 'log'."""
        return self._ax.get_xscale()
    @xscale.setter
    def xscale(self, val):
        if (val != 'linear') and (val != 'log'):
            raise ValueError("Scale must be either 'linear' or 'log'")
        self._ax.set_xscale(val)

    @property
    def ylim(self):
        """(float, float): The plotting range of the y axis. """
        return self._ax.get_ylim()
    @ylim.setter
    def ylim(self, val):
        self._ax.set_ylim(val)

    @property
    def yscale(self):
        """(str): The scale of the y axis, either 'linear' or 'log'. """
        return self._ax.get_yscale()
    @yscale.setter
    def yscale(self, val):
        if (val != 'linear') and (val != 'log'):
            raise ValueError("Scale must be either 'linear' or 'log'")
        self._ax.set_yscale(val)

    def _format_coordinates(self):
        # Set the format of the coordinate readout
        self._ax.format_coord = lambda x, y: ""

        # Set the x minor tick frequency
        minorLocator = AutoMinorLocator()
        self._ax.xaxis.set_minor_locator(minorLocator)

        # Set the y minor tick frequency
        minorLocator = AutoMinorLocator()
        self._ax.yaxis.set_minor_locator(minorLocator)


class GdtCmap():
    """A specialized colormap class that implements an alpha gradient.
    A callback function can also be defined so that when any of the parameters
    are changed, the callback will be executed.
    
    Parameters:
        name (str): Any standard Matplotlib colormap name
        alpha_min (float, optional): The minimum of the alpha gradient, between
                                     (0, 1) inclusive. Default is 0.0.
        alpha_max (float, optional): The maximum of the alpha gradient, between
                                     (0, 1) inclusive. Default is 1.0.
        alpha_scale (str, optional): Set to 'linear' for a linear gradient, 
                                     or to 'log' for a logarithmic gradient.
                                     Default is 'linear'
        callback (function, optional): The callback function to execute when
                                       a parameter is changed.
    """
    def __init__(self, name, alpha_min=0.0, alpha_max=1.0, alpha_scale='linear',
                 callback=None):
        self._name = None
        self._cmap = None
        self._alpha_min = 0.0
        self._alpha_max = 1.0
        self._alpha_scale = 'linear'
        self._callback = None
        
        self.name = name
        if alpha_min > alpha_max:
            raise ValueError('alpha_min must be <= alpha_max')
        self.alpha_min = alpha_min
        self.alpha_max = alpha_max
        self.alpha_scale = alpha_scale
        self.set_callback(callback)
        
    @property
    def alpha_max(self):
        """(float): The maximum of the alpha gradient"""
        return self._alpha_max
    @alpha_max.setter
    def alpha_max(self, val):
        if val < 0.0 or val > 1.0:
            raise ValueError('alpha_max must be >= 0 and <= 1')
        if val < self.alpha_min:
            raise ValueError('alpha_max must be >= alpha_min')
        self._alpha_max = val
        self._update()
        if self._callback is not None:
            self._callback()

    @property
    def alpha_min(self):
        """(float): The minimum of the alpha gradient"""
        if (self.alpha_scale == 'log') and (self._alpha_min == 0):
            return 1e-10
        else:
            return self._alpha_min
    @alpha_min.setter
    def alpha_min(self, val):
        if val < 0.0 or val > 1.0:
            raise ValueError('alpha_min must be >= 0 and <= 1')
        if val > self.alpha_max:
            raise ValueError('alpha_min must be <= alpha_max')
        self._alpha_min = val
        self._update()
        if self._callback is not None:
            self._callback()

    @property
    def alpha_scale(self):
        """(str): The scale of the alpha gradient: linear or log"""
        return self._alpha_scale
    @alpha_scale.setter
    def alpha_scale(self, val):
        if val not in ['linear', 'log']:
            raise ValueError("alpha_scale can only be 'linear' or 'log'")
        self._alpha_scale = val
        self._update()
        if self._callback is not None:
            self._callback()
     
    @property
    def cmap(self):
        """(LinearSegmentedColormap): The custom Matplotlib colormap"""
        return self._cmap
    
    @property
    def name(self):
        """(str): The name of the base colormap"""
        return self._name
    @name.setter
    def name(self, val):
        self._cmap = copy.copy(plt.cm.get_cmap(val))
        self._name = val
        self._update()
        if self._callback is not None:
            self._callback()
    
    def set_callback(self, callback):
        """Set the callback function
        
        Args:
            callback (function): The callback function
        """
        if callback is not None:
            try:
                callback()
            except:
                raise RuntimeError('Not a valid callback function')
        self._callback = callback
    
    def _update(self):
        self._cmap._init()
        self._cmap.set_bad(alpha=0.0)
        
        if self.alpha_scale == 'linear':
            func = np.linspace
        elif self.alpha_scale == 'log':
            func = np.geomspace
        
        self._cmap._lut[:-3,-1] = func(self.alpha_min, self.alpha_max, 
                                       self._cmap.N)
    
    def __repr__(self):
        spaces = ' '*10
        s = '<GdtCmap: {};\n'.format(self.name)
        s += '{0}alpha_min={1:.2f};\n'.format(spaces, self.alpha_min)
        s += '{0}alpha_max={1:.2f};\n'.format(spaces, self.alpha_max)
        s += '{0}alpha_scale={1}>'.format(spaces, self.alpha_scale)
        return s
    

class PlotElement():
    """A base class representing a plot element.  A plot element can be a 
    complex collection of more primitive matplotlib plot elements, but are 
    treated as a single element.
    
    Note:
        This class is not intended to be instantiated directly by the user, 
        rather it is inherited by child plot element objects.

    Parameters:
        alpha (float): The alpha opacity value, between 0 and 1.
        color (str): The color of the plot element
    """
    def __init__(self, color=None, alpha=None):
        self._artists = []
        self._visible = True
        self._color = color
        self._alpha = alpha
        self._kwargs = None

    def __del__(self):
        self.remove()

    @property
    def alpha(self):
        """(float): The alpha opacity value, between 0 and 1."""
        return self._alpha
    @alpha.setter
    def alpha(self, alpha):
        [artist.set_alpha(alpha) for artist in self._artists]
        self._alpha = alpha

    @property
    def artists(self):
        """(list): The object references to the individual matplotlib
        elements """
        return self._artists

    @property
    def color(self):
        """(str): The color of the plot element"""
        return self._color
    @color.setter
    def color(self, color):
        [artist.set_color(color) for artist in self._artists]
        self._color = color

    @property
    def visible(self):
        """(bool): True if the element is shown on the plot, False otherwise"""
        return self._visible

    @property
    def zorder(self):
        """(int): The plot element zorder"""
        return self._artists[0].get_zorder()
    @zorder.setter
    def zorder(self, val):
        for artist in self.artists:
            artist.set_zorder(val)

    def hide(self):
        """Hide the plot element"""
        self._change_visibility(False)

    def remove(self):
        """Remove the plot element"""
        for artist in self._artists:
            try:
                artist.remove()
            except:
                pass

    def show(self):
        """Show the plot element"""
        self._change_visibility(True)

    def toggle(self):
        """Toggle the visibility of the plot element"""
        self._change_visibility(not self._visible)

    def _change_visibility(self, visible):
        for artist in self._artists:
            try:
                artist.set_visible(visible)
            except:
                self._set_visible(artist, visible)
        self._visible = visible

    def _sanitize_artists(self, old_artists):
        """Each artist collection is a collection of artists, and each artist
        is a collection of elements.  Matplotlib isn't exactly consistent
        on how each of these artists are organized for different artist 
        classes, so we have to do some sanitizing
        """
        artists = []
        for artist in old_artists:
            if artist is None:
                continue
            elif isinstance(artist, (tuple, list)):
                if len(artist) == 0:
                    continue
                artists.extend(self._sanitize_artists(artist))
            else:
                artists.append(artist)
        return artists

    def _set_visible(self, artist, visible):
        """In some cases, the artist is a collection of sub-artists and the 
        set_visible() function has not been exposed to the artists, so 
        we must iterate.
        """
        try:
            for subartist in artist.collections:
                subartist.set_visible(visible)
        except:
            pass


class Histo(PlotElement):
    """Plot a rate histogram of either lightcurves or count spectra.
        
    Parameters:
        bins (:class:`~gdt.core.data_primitives.TimeBins` or \
              :class:`~gdt.core.data_primitives.EnergyBins`): 
            The lightcurve or count spectrum histograms
        ax (:class:`matplotlib.axes`): The axis on which to plot
        color (str, optional): The color of the rate histograms
        alpha (float, optional): The alpha of the fill. Default is 1
        label (sr, optional): Label for the histogram
        **kwargs: Other plotting options
    """
    def __init__(self, bins, ax, color='C0', alpha=1.0, label=None, **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(bins, ax)
        self._artists = self._sanitize_artists(artists)
        if label is not None:
            self._artists[0].set_label(label)
    
    @property
    def linestyle(self):
        """(str): The linestyle of the histogram. Default is '-'"""
        return self._artists[0].get_linestyle()
    @linestyle.setter
    def linestyle(self, val):
        for artist in self._artists:
            artist.set_linestyle(val)
            
    @property
    def linewidth(self):
        """(int): The line width of the histogram. Default is 1.5"""
        return self._artists[0].get_linewidth()
    @linewidth.setter
    def linewidth(self, val):
        for artist in self._artists:
            artist.set_linewidth(val)
    
    def _create(self, bins, ax):
        return histo(bins, ax, color=self._color, alpha=self._alpha,
                     **self._kwargs)

    def __repr__(self):
        spaces = ' '*8
        s = '<Histo: color={};\n'.format(self.color) 
        s += "{0}alpha={1};\n{0}linestyle='{2}';\n".format(spaces, self.alpha,
                                                        self.linestyle)
        s += '{0}linewidth={1}>'.format(spaces, self.linewidth)
        return s


class HistoErrorbars(PlotElement):
    """Plot errorbars for lightcurves or count spectra.

    Parameters:
        bins (:class:`~gdt.core.data_primitives.TimeBins` or \
              :class:`~gdt.core.data_primitives.EnergyBins`):
            The lightcurve or count spectrum histograms
        ax (:class:`matplotlib.axes`): The axis on which to plot
        color (str, optional): The color of the rate histograms
        alpha (float, optional): The alpha of the fill. Default is 1
        **kwargs: Other plotting options
    """
    def __init__(self, bins, ax, color=None, alpha=None, **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(bins, ax)
        self._artists = self._sanitize_artists(artists)

    @property
    def linewidth(self):
        """(int): The line width of the errorbars. Default is 1.5"""
        return self._artists[0].get_linewidth()[0]
    @linewidth.setter
    def linewidth(self, val):
        for artist in self._artists:
            artist.set_linewidth(val)
    
    def _create(self, bins, ax):
        return histo_errorbars(bins, ax, color=self._color, alpha=self._alpha,
                               **self._kwargs)

    def __repr__(self):
        spaces = ' '*17
        s = '<HistoErrorbars: color={0};\n{1}'.format(self.color, spaces)
        s += 'alpha={0};\n{1}linewidth={2}>'.format(self.alpha, spaces,
                                                    self.linewidth)
        return s


class HistoFilled(PlotElement):
    """Plot a filled histogram.

    Parameters:
        bins (:class:`~gdt.core.data_primitives.TimeBins` or \
              :class:`~gdt.core.data_primitives.EnergyBins`):
            The lightcurve or count spectrum histograms
        ax (:class:`matplotlib.axes`): The axis on which to plot
        color (str, optional): The color of the rate histograms
        alpha (float, optional): The alpha of the edges. Default is 1
        fill_alpha (float, optional): The alpha of the fill. Default is 0.2
        **kwargs: Other plotting options
    """
    def __init__(self, bins, ax, color=None, alpha=None, fill_alpha=0.2, 
                 **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(bins, ax)
        self._artists = self._sanitize_artists(artists)
        self.fill_alpha = fill_alpha

    @property
    def alpha(self):
        """(float): The alpha opacity value of the edge"""
        return self._artists[0].get_alpha()
    @alpha.setter
    def alpha(self, val):
        for artist in self._artists[:-1]:
            artist.set_alpha(val)

    @property
    def fill_alpha(self):
        """(float): The alpha opacity value of the fill"""
        return self._artists[-1].get_alpha()
    @fill_alpha.setter
    def fill_alpha(self, val):
        self._artists[-1].set_alpha(val)

    @property
    def linestyle(self):
        """(str): The linestyle of the histogram. Default is '-'"""
        return self._artists[0].get_linestyle()
    @linestyle.setter
    def linestyle(self, val):
        for artist in self._artists:
            artist.set_linestyle(val)
            
    @property
    def linewidth(self):
        """(int): The line width of the histogram. Default is 1.5"""
        return self._artists[0].get_linewidth()
    @linewidth.setter
    def linewidth(self, val):
        for artist in self._artists:
            artist.set_linewidth(val)
    
    def _create(self, bins, ax):
        return histo_filled(bins, ax, color=self._color,
                            fill_alpha=self._alpha,
                            **self._kwargs)

    def __repr__(self):
        spaces = ' '*14
        s = '<HistoFilled: color={0};\n{1}'.format(self.color, spaces) 
        s += 'alpha={0};\n{1}fill_alpha={2};\n{1}'.format(self.alpha, spaces,
                                                          self.fill_alpha)
        s += "linestyle='{0}';\n{1}linewidth={2}>".format(self.linestyle, 
                                                          spaces,self.linewidth)
        return s


class LightcurveBackground(PlotElement):
    """Plot a lightcurve background model with an error band.

    Parameters:
        backrates (:class:`~gdt.background.primitives.BackgroundRates`):
            The background rates object integrated over energy.  If there is 
            more than one remaining energy channel, the background will be 
            integrated over the remaining energy channels.
        ax (:class:`matplotlib.axes`): The axis on which to plot
        color (str, optional): The color of the background.
        alpha (float, optional): The alpha of the background central line.
        band_alpha (float, optional): The alpha of the background uncertainty. 
        **kwargs: Other plotting options
    """
    def __init__(self, backrates, ax, color=None, alpha=None, band_alpha=None, 
                 **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(backrates, ax)
        self._artists = self._sanitize_artists(artists)
        self.band_alpha = band_alpha

    @property
    def alpha(self):
        """(float): The alpha opacity value"""
        return self._artists[0].get_alpha()
    @alpha.setter
    def alpha(self, alpha):
        self._artists[0].set_alpha(alpha)

    @property
    def band_alpha(self):
        """(float): The opacity of the uncertainty band"""
        return self._artists[1].get_alpha()
    @band_alpha.setter
    def band_alpha(self, alpha):
        return self._artists[1].set_alpha(alpha)

    @property
    def linestyle(self):
        """(str): The linestyle of the background Default is '-'"""
        return self._artists[0].get_linestyle()
    @linestyle.setter
    def linestyle(self, val):
        self._artists[0].set_linestyle(val)
            
    @property
    def linewidth(self):
        """(int): The line width of the background. Default is 1.5"""
        return self._artists[0].get_linewidth()
    @linewidth.setter
    def linewidth(self, val):
        self._artists[0].set_linewidth(val)
    
    def _create(self, backrates, ax):
        return lightcurve_background(backrates, ax, cent_color=self.color,
                                     err_color=self.color, 
                                     cent_alpha=self._alpha, 
                                     err_alpha=self._alpha,
                                     **self._kwargs)

    def __repr__(self):
        spaces = ' '*23
        s = '<LightcurveBackground: color={0};\n{1}'.format(self.color, spaces) 
        s += 'alpha={0};\n{1}band_alpha={2};\n{1}'.format(self.alpha, spaces,
                                                          self.band_alpha)
        s += "linestyle='{0}';\n{1}linewidth={2}>".format(self.linestyle, 
                                                          spaces,self.linewidth)
        return s


class SpectrumBackground(PlotElement):
    """Plot a count spectrum background model with an error band.

    Parameters:
        backspec (:class:`~gdt.background.primitives.BackgroundSpectrum`):
            The background sepctrum object integrated over time
        ax (:class:`matplotlib.axes`): The axis on which to plot
        cent_alpha (float, optional): The alpha of the background centroid line. 
                                      Default is 1
        cent_color (str, optional): The color of the background centroid line
        err_alpha (float, optional): The alpha of the background uncertainty. 
                                     Default is 1
        err_color (str, optional): The color of the background uncertainty
        
        color (str, optional): The color of the background. If set, overrides 
                               ``cent_color`` and ``err_color``.
        alpha (float, optional): The alpha of the background. If set, overrides 
                                 ``cent_alpha`` and ``err_alpha``
        **kwargs: Other plotting options
    """
    def __init__(self, backspec, ax, color=None, alpha=None, band_alpha=None, 
                **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(backspec, ax)
        self._artists = self._sanitize_artists(artists)
        self.band_alpha = band_alpha

    @property
    def alpha(self):
        """(float): The alpha opacity value"""
        return self._artists[0].get_alpha()
    @alpha.setter
    def alpha(self, alpha):
        self._artists[0].set_alpha(alpha)

    @property
    def band_alpha(self):
        """(float): The opacity of the uncertainty band"""
        return self._artists[1].get_alpha()
    @band_alpha.setter
    def band_alpha(self, alpha):
        return self._artists[1].set_alpha(alpha)

    @property
    def linestyle(self):
        """(str): The linestyle of the background Default is '-'"""
        return self._artists[0].get_linestyle()
    @linestyle.setter
    def linestyle(self, val):
        self._artists[0].set_linestyle(val)
            
    @property
    def linewidth(self):
        """(int): The line width of the background. Default is 1.5"""
        return self._artists[0].get_linewidth()
    @linewidth.setter
    def linewidth(self, val):
        self._artists[0].set_linewidth(val)
    
    def _create(self, backspec, ax):
        return spectrum_background(backspec, ax, cent_color=self.color,
                                   err_color=self.color, cent_alpha=self._alpha,
                                   err_alpha=self._alpha, **self._kwargs)

    def __repr__(self):
        spaces = ' '*21
        s = '<SpectrumBackground: color={0};\n{1}'.format(self.color, spaces) 
        s += 'alpha={0};\n{1}band_alpha={2};\n{1}'.format(self.alpha, spaces,
                                                          self.band_alpha)
        s += "linestyle='{0}';\n{1}linewidth={2}>".format(self.linestyle, 
                                                          spaces,self.linewidth)
        return s


class Heatmap(PlotElement):
    """Plot a general heatmap (color gradient). By default, the heatmap opacity 
    will scale from 0 (fully transparent) to 1 (fully opaque) corresponding to 
    the colormap value, creating an opacity gradient. This behavior can be 
    adjusted by setting ``alpha_min`` and ``alpha_max``.
    
    Parameters:
        x (np.array): The x-coordinate array of the heatmap grid
        y (np.array): The y-coordinate array of the heatmap grid
        heatmap (np.array): The heatmap array, of shape (``x.size``, ``y.size``)
        ax (:class:`matplotlib.axes`): The axis on which to plot
        colorbar (bool, optional): If True, add the colorbar scale. 
                                   Default is True.
        color (:class:`GdtCmap`): The colormap of the heatmap. 
                                  Default is 'viridis'
        norm (:class:`matplotlib.colors.Normalize` or similar, optional):
            The normalization used to scale the colormapping to the heatmap 
            values. This can be initialized by the defined matplotlib 
            normalizations or a custom normalization. 
        num_contours (int, optional): The number of contour lines to draw.
                                      Default is 100.
        **kwargs: Other plotting options
    """
    def __init__(self, x, y, heatmap, ax, colorbar=True, color=GdtCmap('viridis'),
                 alpha=None, norm=None, num_contours=100, **kwargs):
        
        # store data for replotting if needed
        self._x = x
        self._y = y
        self._heatmap = heatmap
        self._ax = ax
        self._colorbar = colorbar
        self._num_contours = num_contours
        
        # set the normalization, and ensure that it is scaled to the data range
        if norm is None:
            norm = PowerNorm(gamma=0.3)
        self._norm = norm
        self._norm.vmin = self._heatmap.min()
        self._norm.vmax = self._heatmap.max()

        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create()
        self._artists = self._sanitize_artists(artists)
        
        # set the colormap
        self.color = color
        self.color.set_callback(self._artists[0].changed)
    
    @property
    def color(self):
        """(GdtCmap): The colormap"""
        return self._color
    @color.setter
    def color(self, color):
        if not isinstance(color, GdtCmap):
            raise TypeError('color must be of type GdtCmap')
        self._color = color
        self._artists[0].set_cmap(color.cmap)

    @property
    def colorbar(self):
        """(matplotlib.colorbar.Colorbar): The colorbar object"""
        if self._colorbar:
            return self._artists[-1]
    
    @property
    def norm(self):
        """(matplotlib.colors.Normalize or similar): The colormap normalization"""
        return self._norm
    @norm.setter
    def norm(self, norm):
        # mark FIXME: This creates a new colorbar.  Need to force it to update
        # or go through _make_colorbar()
        norm.vmin = self._heatmap.min()
        norm.vmax = self._heatmap.max()
        self._artists[0].set_norm(norm)
        self._norm = norm

    @property
    def num_contours(self):
        """(int): Number of plot contours"""
        # mark TODO: Add a setter, which will have to do a replot
        return self._num_contours
    
    def _create(self):
        refs = response_matrix(self._x, self._y, self._heatmap, self._ax, 
                               cmap=self._color.cmap, norm=self._norm, 
                               num_contours=self._num_contours, **self._kwargs)
        refs = [refs]
        if self._colorbar:
            refs.append(self._make_colorbar(self._ax, refs[0]))

        return refs

    def _make_colorbar(self, ax, artist):
        # scientific format for the colorbar
        def sci_fmt(x, pos):
            if x >= 1e-3:
                return r'{:.3f}'.format(x)
            else:
                a, b = '{:.0e}'.format(x).split('e')
                b = int(b)
                return r'${} \times 10^{{{}}}$'.format(a, b)
        
        # vertical colorbar
        cb = plt.colorbar(artist, ax=ax, aspect=10.0, pad=0.01, panchor=False,
                          cmap=self._color.cmap, norm=self._norm,
                          orientation='vertical',
                          format=ticker.FuncFormatter(sci_fmt))
        cb.ax.tick_params(labelsize=10)
        cb.set_label(r'Effective Area (cm$^2$)', fontsize=12)
        cb.patch.set_facecolor('black')
        cb.draw_all()
        return cb

    def __repr__(self):
        spaces = ' ' * 10
        s = "<Heatmap: color='{0}';\n{1}".format(self.color.name, spaces) 
        s += 'norm={0};\n{1}'.format(self.norm.__class__.__name__, spaces)
        s += 'num_contours={0};\n{1}'.format(self.num_contours, spaces)
        s += 'colorbar={0}>'.format(self._colorbar)
        return s


class EffectiveArea(PlotElement):
    """Plot a histogram of the effective area of a detector response.
    
    Parameters:
        bins (:class:`~gdt.core.data_primitives.Bins`): The effective area 
                                                        histograms
        ax (:class:`matplotlib.axes`): The axis on which to plot
        color (str, optional): The color of the rate histograms
        alpha (float, optional): The alpha of the histogram. Default is 1
        **kwargs: Other plotting options
    """
    def __init__(self, bins, ax, color='C0', alpha=1.0, orientation='vertical',
                 **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        self._orientation = orientation
        artists = self._create(bins, ax)
        self._artists = self._sanitize_artists(artists)

    @property
    def linestyle(self):
        """(str): The linestyle of the histogram. Default is '-'"""
        return self._artists[0].get_linestyle()
    @linestyle.setter
    def linestyle(self, val):
        for artist in self._artists:
            artist.set_linestyle(val)
            
    @property
    def linewidth(self):
        """(int): The line width of the histogram. Default is 1.5"""
        return self._artists[0].get_linewidth()
    @linewidth.setter
    def linewidth(self, val):
        for artist in self._artists:
            artist.set_linewidth(val)
    
    def _create(self, bins, ax):
        return effective_area(bins, ax, color=self._color, alpha=self._alpha,
                              orientation=self._orientation, **self._kwargs)

    def __repr__(self):
        spaces = ' ' * 16
        s = '<EffectiveArea: color={};\n'.format(self.color) 
        s += "{0}alpha={1};\n{0}linestyle='{2}';\n".format(spaces, self.alpha,
                                                        self.linestyle)
        s += '{0}linewidth={1}>'.format(spaces, self.linewidth)
        return s


class SAA(PlotElement):
    """Plot the SAA polygon on the Earth.

    Parameters:
        saa_def (:class:`~gdt.core.geomagnetic.SouthAtlanticAnomaly`):
            Object containing the SAA boundary definition
        proj (Cartopy Projection): The Cartopy projection
        ax (:class:`matplotlib.axes`): The axis on which to plot
        color (str, optional): The color of the SAA polygon
        alpha (float, optional): The alpha of the interior of the polygon
        **kwargs: Other plotting options
    """
    def __init__(self, saa_def, proj, ax, color='darkred', alpha=0.4, **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(saa_def, proj, ax)
        self._artists = self._sanitize_artists(artists)
        self._path = artists[0].get_path()

    @property
    def fill(self):
        """(bool): True if the polygon is filled, False otherwise"""
        return [artist.get_fill() for artist in self.artists \
                if artist.__class__.__name__ == 'Polygon'][0]
    @fill.setter
    def fill(self, val):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Polygon':
                artist.set_fill(val)

    @property
    def hatch(self):
        """(str): The hatching string"""
        return [artist.get_hatch() for artist in self.artists \
                if artist.__class__.__name__ == 'Polygon'][0]
    @hatch.setter
    def hatch(self, val):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Polygon':
                artist.set_hatch(val)

    @property
    def linestyle(self):
        """(str): The linestyle of the edge"""
        return [artist.get_linestyle() for artist in self.artists \
                if artist.__class__.__name__ == 'Polygon'][0]
    @linestyle.setter
    def linestyle(self, val):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Polygon':
                artist.set_linestyle(val)

    @property
    def linewidth(self):
        """(int): The line width of the edge"""
        return [artist.get_linewidth() for artist in self.artists \
                if artist.__class__.__name__ == 'Polygon'][0]
    @linewidth.setter
    def linewidth(self, val):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Polygon':
                artist.set_linewidth(val)

    def in_saa(self, lon, lat):
        """Check if a point or points is inside the SAA
    
        Args:
            lon (np.array): longitudes
            lat (np.array): latitudes
    
        Returns:
            np.array: Boolean array for each point where True indicates the \
                      point is in the SAA.
        """
        pts = np.array((lon.ravel(), lat.ravel())).T
        mask = self._path.contains_points(pts)
        return mask

    def _create(self, saa_def, proj, ax):
        artists = saa_polygon(saa_def.latitude, saa_def.longitude, proj, 
                              color=self._color, alpha=self._alpha,
                              **self._kwargs)
        ax.add_patch(artists)
        return [artists]

    def __repr__(self):
        spaces = ' '*6
        s = "<SAA: color='{}';\n".format(self.color) 
        s += "{0}alpha={1};\n".format(spaces, self.alpha)
        s += "{0}linestyle='{1}';\n".format(spaces, self.linestyle)
        s += '{0}linewidth={1}>'.format(spaces, self.linewidth)
        return s


class EarthLine(PlotElement):
    """Plot a line on the Earth.
    
    Parameters:
        lat (np.array): The latitude values of the line
        lon (np.array): The longitude values of the line
        proj (GeoAxesSubplot): The Cartopy projection
        color (str, optional): The color of the line
        alpha (float, optional): The alpha opacity of the line
        **kwargs: Other plotting options
    """
    def __init__(self, lat, lon, proj, color='black', alpha=0.4, **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(lat, lon, proj)
        self._artists = self._sanitize_artists(artists)

    @property
    def linestyle(self):
        """(str): The linestyle of the histogram. Default is '-'"""
        return self._artists[0].get_linestyle()
    @linestyle.setter
    def linestyle(self, val):
        for artist in self._artists:
            artist.set_linestyle(val)
            
    @property
    def linewidth(self):
        """(int): The line width of the histogram. Default is 1.5"""
        return self._artists[0].get_linewidth()
    @linewidth.setter
    def linewidth(self, val):
        for artist in self._artists:
            artist.set_linewidth(val)

    def _create(self, lat, lon, proj):
        artists = earth_line(lat, lon, proj, color=self._color,
                             alpha=self._alpha, **self._kwargs)
        return artists

    def __repr__(self):
        spaces = ' ' * 12
        s = '<EarthLine: color={};\n'.format(self.color) 
        s += "{0}alpha={1};\n{0}linestyle='{2}';\n".format(spaces, self.alpha,
                                                           self.linestyle)
        s += '{0}linewidth={1}>'.format(spaces, self.linewidth)
        return s


class EarthPoints(PlotElement):
    """Plot a point or set of points on the Earth.
    
    Parameters:
        lat (np.array): The latitude values
        lon (np.array): The longitude values
        proj (GeoAxesSubplot): The Cartopy projection
        color (str, optional): The color of the line
        alpha (float, optional): The alpha opacity of the line
        sizes (list, optional): The sizes of the plot points. Default is 10.
        **kwargs: Other plotting options
    """
    def __init__(self, lat, lon, proj, color='black', alpha=1.0, sizes=10, 
                 **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(lat, lon, proj)
        self._artists = self._sanitize_artists(artists)
        self._lat = lat
        self._lon = lon
        self.sizes = sizes

    @property
    def coordinates(self):
        """(list of str): The formatted coordinate list of the points"""
        coords = self._to_geocoords()
        return coords

    @property
    def num_points(self):
        """(int): Number of plotted points"""
        return len(self.sizes)
    
    @property
    def sizes(self):
        """(list): The size of each plotted point"""
        try:
            return self._artists[0].get_sizes()
        except:
            return [1]
    @sizes.setter
    def sizes(self, val):
        if isinstance(val, (int, float)):
            val = [val] * self.num_points
        try:
            self._artists[0].set_sizes(val)
        except:
            pass

    def _create(self, lat, lon, proj):
        artists = earth_points(lat, lon, proj, color=self._color,
                               alpha=self._alpha, **self._kwargs)
        return artists

    def _to_geocoords(self):
        lat_formatter = lambda x: "%.2fN" % (x) if x > 0 else "%.2fS" % (-x)
        lon_formatter = lambda x: "%.2fE" % (x) if x > 0 else "%.2fW" % (-x)
        lat = lat_formatter(self._lat)
        lon = lon_formatter(self._lon)
        return (lat, lon)

    def __repr__(self):
        spaces = ' ' * 14
        s = '<EarthPoints: {0} points;\n{1}'.format(self.num_points, spaces)
        s += "color='{0}';\n{1}alpha={2}>".format(self.color, spaces, self.alpha) 
        return s


class SkyPoints(PlotElement):
    """Plot a set of points on the sky in equatorial, galactic, or spacecraft 
    coordinates.
    
    Parameters:
        x (float or np.array): 
            The azimuthal coordinate, in degrees
        y (float or np.array):
            The polar coordinate, in degrees
        ax (:class:`matplotlib.axes`): The axis on which to plot
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', 'spacecraft'.
                               Default is 'equatorial'
        color (str, optional): The color of the points
        alpha (float, optional): The alpha opacity of the points
        **kwargs: Other plotting options
    """
    def __init__(self, x, y, ax, flipped=True, frame='equatorial', color='C0',
                 alpha=1.0, **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(x, y, ax, flipped, frame)
        self._artists = self._sanitize_artists(artists)

    @property
    def num_points(self):
        """(int): The number of points plotted"""
        return len(self.sizes)
    
    @property
    def sizes(self):
        """(list): The size of each plot point"""
        return self._artists[0].get_sizes()
    @sizes.setter
    def sizes(self, val):
        if isinstance(val, (int, float)):
            val = [val] * self.num_points
        self._artists[0].set_sizes(val)
            
    def _create(self, x, y, ax, flipped, frame):
        ref = sky_point(x, y, ax, flipped=flipped, frame=frame,
                        color=self._color, alpha=self._alpha, **self._kwargs)
        return [ref]

    def __repr__(self):
        spaces = ' '*12
        s = '<SkyPoints: {0} points;\n{1}'.format(self.num_points, spaces)
        s += "color='{0}';\n{1}alpha={2}>".format(self.color, spaces, self.alpha) 
        return s


class SkyLine(PlotElement):
    """Plot a line on the sky in either equatorial, galactic, or spacecraft 
    coordinates.

    Parameters:
        x (np.array): The azimuthal coordinate, in degrees
        y (np.array): The polar coordinate, in degrees
        ax (:class:`matplotlib.axes`): The axis on which to plot
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', 'spacecraft'.
                               Default is 'equatorial'
        color (str, optional): The color of the line
        alpha (float, optional): The alpha opacity of the line
        **kwargs: Other plotting options
    """
    def __init__(self, x, y, ax, flipped=True, frame='equatorial', color='C0',
                 alpha=1.0, **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(x, y, ax, flipped, frame)
        self._artists = self._sanitize_artists(artists)

    @property
    def linestyle(self):
        """(str): The linestyle of the histogram. Default is '-'"""
        return self._artists[0].get_linestyle()
    @linestyle.setter
    def linestyle(self, val):
        for artist in self._artists:
            artist.set_linestyle(val)
            
    @property
    def linewidth(self):
        """(int): The line width of the histogram. Default is 1.5"""
        return self._artists[0].get_linewidth()
    @linewidth.setter
    def linewidth(self, val):
        for artist in self._artists:
            artist.set_linewidth(val)

    def _create(self, x, y, ax, flipped, frame):
        refs = sky_line(x, y, ax, flipped=flipped, frame=frame,
                        color=self._color, alpha=self._alpha, **self._kwargs)
        return refs

    def __repr__(self):
        spaces = ' '*10
        s = '<SkyLine: color={};\n'.format(self.color) 
        s += "{0}alpha={1};\n{0}linestyle='{2}';\n".format(spaces, self.alpha,
                                                        self.linestyle)
        s += '{0}linewidth={1}>'.format(spaces, self.linewidth)
        return s


class SkyCircle(PlotElement):
    """Plot a circle on the sky in equatorial, galactic, or spacecraft 
    coordinates.

    Parameters:
        x (float): The azimuthal coordinate of the center, in degrees
        y (float): The polar coordinate of the center, in degrees
        radius (float): The radius of the circle, in degrees
        ax (:class:`matplotlib.axes`): The axis on which to plot
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', 'spacecraft'.
                               Default is 'equatorial'
        face_color (str, optional): The color of the circle fill
        face_alpha (float, optional): The alpha opacity of the circle fill
        edge_color (str, optional): The color of the circle edge
        edge_alpha (float, optional): The alpha opacity of the circle edge
        color (str, optional): The color of the circle. If set, overrides 
                               ``face_color`` and ``edge_color``.
        alpha (float, optional): The alpha of the circle. If set, overrides 
                               ``face_alpha`` and ``edge_alpha``.
        **kwargs: Other plotting options
    """
    def __init__(self, x, y, radius, ax, flipped=True, frame='equatorial', 
                 color='C0', alpha=1.0, face_color=None, face_alpha=None, 
                 edge_color=None, edge_alpha=None, **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        # color and alpha act as setting the global color values for the object
        if color is not None:
            face_color = color
            edge_color = color
        if alpha is not None:
            face_alpha = alpha
            edge_alpha = alpha
        self._face_alpha = face_alpha
        self._edge_alpha = edge_alpha
        self._face_color = face_color
        self._edge_color = edge_color

        artists = self._create(x, y, radius, ax, flipped, frame)
        self._artists = self._sanitize_artists(artists)

    @property
    def alpha(self):
        """(float): The alpha opacity value"""
        return self._alpha
    @alpha.setter
    def alpha(self, alpha):
        self.face_alpha = alpha
        self.edge_alpha = alpha
        self._alpha = alpha

    @property
    def color(self):
        """(str): The color of the circle"""
        return self._color
    @color.setter
    def color(self, color):
        self.face_color = color
        self.edge_color = color
        self._color = color

    @property
    def edge_alpha(self):
        """(float): The alpha opacity of the circle edge"""
        return self._edge_alpha
    @edge_alpha.setter
    def edge_alpha(self, alpha):
        edge = colorConverter.to_rgba(self._edge_color, alpha=alpha)
        [artist.set_edgecolor(edge) for artist in self._artists \
         if hasattr(artist, 'set_edgecolor')]
        self._edge_alpha = alpha

    @property
    def edge_color(self):
        """(str): The color of the circle edge"""
        return self._edge_color
    @edge_color.setter
    def edge_color(self, color):
        edge = colorConverter.to_rgba(color, alpha=self._edge_alpha)
        [artist.set_edgecolor(edge) for artist in self._artists \
         if hasattr(artist, 'set_edgecolor')]
        self._edge_color = color

    @property
    def face_alpha(self):
        """(float): The alpha opacity of the circle fill"""
        return self._face_alpha
    @face_alpha.setter
    def face_alpha(self, alpha):
        face = colorConverter.to_rgba(self._face_color, alpha=alpha)
        [artist.set_facecolor(face) for artist in self._artists \
         if hasattr(artist, 'set_facecolor')]
        self._face_alpha = alpha

    @property
    def face_color(self):
        """(str): The color of the circle fill"""
        return self._face_color
    @face_color.setter
    def face_color(self, color):
        face = colorConverter.to_rgba(color, alpha=self._face_alpha)
        [artist.set_facecolor(face) for artist in self._artists \
         if hasattr(artist, 'set_facecolor')]
        self._face_color = color

    @property
    def fill(self):
        """(bool): True if the circle is filled, False otherwise"""
        return self._artists[0].get_fill()
    @fill.setter
    def fill(self, val):
        for artist in self._artists:
            try:
                artist.set_fill(val)
            except:
                pass
        
    @property
    def hatch(self):
        """(str): The hatching string"""
        return self._artists[0].get_hatch()
    @hatch.setter
    def hatch(self, val):
        for artist in self._artists:
            try:
                artist.set_hatch(val)
            except:
                pass

    @property
    def linestyle(self):
        """(str): The linestyle of the edge"""
        return self._artists[0].get_linestyle()
    @linestyle.setter
    def linestyle(self, val):
        for artist in self._artists:
            try:
                artist.set_linestyle(val)
            except:
                pass
            
    @property
    def linewidth(self):
        """(int): The line width of the edge"""
        return self._artists[0].get_linewidth()
    @linewidth.setter
    def linewidth(self, val):
        for artist in self._artists:
            try:
                artist.set_linewidth(val)
            except:
                pass

    def _create(self, x, y, radius, ax, flipped, frame):
        refs = sky_circle(x, y, radius, ax, flipped=flipped, frame=frame,
                          face_color=self._face_color,
                          face_alpha=self._face_alpha,
                          edge_color=self._edge_color,
                          edge_alpha=self._edge_alpha,
                          **self._kwargs)
        return refs

    def __repr__(self):
        spaces = ' '*12
        s = '<SkyCircle: face_color={};\n'.format(self.face_color) 
        s += "{0}face_alpha={1};\n".format(spaces, self.face_alpha)
        s += "{0}edge_color={1};\n".format(spaces, self.edge_color)
        s += "{0}edge_alpha={1};\n".format(spaces, self.edge_alpha)
        s += "{0}linestyle='{1}';\n".format(spaces, self.linestyle)
        s += '{0}linewidth={1}>'.format(spaces, self.linewidth)
        return s


class SkyPolygon(PlotElement):
    """Plot a polygon on the sky in equatorial, galactic, or spacecraft 
    coordinates.

    Parameters:
        x (float): The azimuthal coordinate, in degrees
        y (float): The polar coordinate, in degrees
        ax (:class:`matplotlib.axes`): The axis on which to plot
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', 'spacecraft'.
                               Default is 'equatorial'
        face_color (str, optional): The color of the polygon fill
        face_alpha (float, optional): The alpha opacity of the polygon fill
        edge_color (str, optional): The color of the polygon edge
        edge_alpha (float, optional): The alpha opacity of the polygon edge
        color (str, optional): The color of the polygon. If set, overrides 
                               ``face_color`` and ``edge_color``.
        alpha (float, optional): The alpha of the polygon. If set, overrides 
                               ``face_alpha`` and ``edge_alpha``.
        **kwargs: Other plotting options
     """
    def __init__(self, x, y, ax, flipped=True, frame='equatorial', color='C0',
                 alpha=1.0, face_color=None, face_alpha=None, edge_color=None,
                 edge_alpha=None, **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        # color and alpha act as setting the global color values for the object
        if color is not None:
            face_color = color
            edge_color = color
        if alpha is not None:
            face_alpha = alpha
            edge_alpha = alpha
        self._face_alpha = face_alpha
        self._edge_alpha = edge_alpha
        self._face_color = face_color
        self._edge_color = edge_color

        artists = self._create(x, y, ax, flipped, frame)
        self._artists = self._sanitize_artists(artists)

    @property
    def alpha(self):
        """(float): The alpha opacity value"""
        return self._alpha
    @alpha.setter
    def alpha(self, alpha):
        self.face_alpha = alpha
        self.edge_alpha = alpha
        self._alpha = alpha

    @property
    def color(self):
        """(str): The color of the polygon"""
        return self._color
    @color.setter
    def color(self, color):
        self.face_color = color
        self.edge_color = color
        self._color = color

    @property
    def edge_alpha(self):
        """(float): The alpha opacity of the polygon edge"""
        return self._edge_alpha
    @edge_alpha.setter
    def edge_alpha(self, alpha):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Line2D':
                artist.set_alpha(alpha)
        self._edge_alpha = alpha

    @property
    def edge_color(self):
        """(str): The color of the polygon edge"""
        return self._edge_color
    @edge_color.setter
    def edge_color(self, color):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Line2D':
                artist.set_color(color)
        self._edge_color = color

    @property
    def face_alpha(self):
        """(float): The alpha opacity of the polygon fill"""
        return self._face_alpha
    @face_alpha.setter
    def face_alpha(self, alpha):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Polygon':
                artist.set_alpha(alpha)
        self._face_alpha = alpha

    @property
    def face_color(self):
        """(str): The color of the polygon fill"""
        return self._face_color
    @face_color.setter
    def face_color(self, color):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Polygon':
                artist.set_color(color)
        self._face_color = color

    @property
    def fill(self):
        """(bool): True if the polygon is filled, False otherwise"""
        return [artist.get_fill() for artist in self.artists \
                if artist.__class__.__name__ == 'Polygon'][0]
    @fill.setter
    def fill(self, val):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Polygon':
                artist.set_fill(val)

    @property
    def hatch(self):
        """(str): The hatching string"""
        return [artist.get_hatch() for artist in self.artists \
                if artist.__class__.__name__ == 'Polygon'][0]
    @hatch.setter
    def hatch(self, val):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Polygon':
                artist.set_hatch(val)

    @property
    def linestyle(self):
        """(str): The linestyle of the edge"""
        return [artist.get_linestyle() for artist in self.artists \
                if artist.__class__.__name__ == 'Line2D'][0]
    @linestyle.setter
    def linestyle(self, val):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Line2D':
                artist.set_linestyle(val)
            
    @property
    def linewidth(self):
        """(int): The line width of the edge"""
        return [artist.get_linewidth() for artist in self.artists \
                if artist.__class__.__name__ == 'Line2D'][0]
    @linewidth.setter
    def linewidth(self, val):
        for artist in self.artists:
            if artist.__class__.__name__ == 'Line2D':
                artist.set_linewidth(val)
        
    def _create(self, x, y, ax, flipped, frame):
        refs = sky_polygon(x, y, ax, flipped=flipped, frame=frame,
                           face_color=self._face_color,
                           face_alpha=self._face_alpha,
                           edge_color=self._edge_color,
                           edge_alpha=self._edge_alpha,
                           **self._kwargs)
        return refs

    def __repr__(self):
        spaces = ' '*13
        s = '<SkyPolygon: face_color={};\n'.format(self.face_color) 
        s += "{0}face_alpha={1};\n".format(spaces, self.face_alpha)
        s += "{0}edge_color={1};\n".format(spaces, self.edge_color)
        s += "{0}edge_alpha={1};\n".format(spaces, self.edge_alpha)
        s += "{0}linestyle='{1}';\n".format(spaces, self.linestyle)
        s += '{0}linewidth={1}>'.format(spaces, self.linewidth)
        return s


class SkyAnnulus(PlotElement):
    """Plot an annulus on the sky in equatorial, galactic, or spacecraft 
    coordinates.
    
    Parameters:
        x (float): The azimuthal center, in degrees
        y (float): The polar center, in degrees
        radius (float): The radius in degrees, defined as the angular distance 
                        from the center to the middle of the width of the 
                        annulus band.
        width (float): The width onf the annulus in degrees
        ax (:class:`matplotlib.axes`): The axis on which to plot
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', 'spacecraft'.
                               Default is 'equatorial'
        color (str, optional): The color of the annulus.
        alpha (float, optional): The opacity of the annulus.
        **kwargs: Other plotting options
    """
    def __init__(self, x, y, radius, width, ax, flipped=True, 
                 frame='equatorial', color='C0', alpha=1.0, **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs

        artists = self._create(x, y, radius, width, ax, flipped, frame)
        self._artists = self._sanitize_artists(artists)

    @property
    def color(self):
        """(str): The color of the annulus"""
        return self._color
    @color.setter
    def color(self, color):
        edge = colorConverter.to_rgba(color, alpha=1.0)
        face = colorConverter.to_rgba(color, alpha=self.alpha)
        [artist.set_edgecolor(edge) for artist in self._artists]
        [artist.set_facecolor(face) for artist in self._artists]
        self._color = color

    @property
    def hatch(self):
        """(str): The hatching string"""
        return self._artists[0].get_hatch()
    @hatch.setter
    def hatch(self, val):
        for artist in self._artists:
            artist.set_hatch(val)

    @property
    def linestyle(self):
        """(str): The linestyle of the edge"""
        return self._artists[0].get_linestyle()
    @linestyle.setter
    def linestyle(self, val):
        for artist in self._artists:
            artist.set_linestyle(val)
            
    @property
    def linewidth(self):
        """(int): The line width of the edge"""
        return self._artists[0].get_linewidth()
    @linewidth.setter
    def linewidth(self, val):
        for artist in self._artists:
            artist.set_linewidth(val)

    def _create(self, x, y, radius, width, ax, flipped, frame):
        refs = sky_annulus(x, y, radius, width, ax, flipped=flipped,
                           frame=frame, color=self._color, alpha=self._alpha,
                           **self._kwargs)
        return refs

    def __repr__(self):
        spaces = ' '*13
        s = '<SkyAnnulus: color={};\n'.format(self.color) 
        s += "{0}alpha={1};\n".format(spaces, self.alpha)
        s += "{0}linestyle='{1}';\n".format(spaces, self.linestyle)
        s += '{0}linewidth={1}>'.format(spaces, self.linewidth)
        return s


class Sun(PlotElement):
    """Plot the sun on the sky in equatorial, galactic, or spacecraft
    coordinates

    Parameters:
        x (float): The azimuthal coordinate, in degrees
        y (float): The polar coordinate, in degrees
        ax (:class:`matplotlib.axes`): The axis on which to plot
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', 'spacecraft'.
                               Default is 'equatorial'
        alpha (float, optional): The opacity of the annulus.
        **kwargs: Other plotting options
    """
    def __init__(self, x, y, ax, flipped=True, frame='equatorial', alpha=1.0,
                 **kwargs):
        super().__init__(color=None, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(x, y, ax, flipped, frame)
        self._artists = self._sanitize_artists(artists)
            
    @property
    def size(self):
        """(float): The size of the plot point"""
        return self.artists[0].get_sizes()[0]
    @size.setter
    def size(self, val):
        frac = val/self.size
        self.artists[0].set_sizes([val])
        face_size = self.artists[1].get_sizes()[0]
        self.artists[1].set_sizes([frac * face_size])
     
    def _create(self, x, y, ax, flipped, frame):
        # :)
        marker1 = "$\u263c$"
        marker2 = "$\u263b$"

        edge = colorConverter.to_rgba('gold', alpha=self._alpha)
        face = colorConverter.to_rgba('yellow', alpha=0.75 * self._alpha)
        point1 = sky_point(x, y, ax, marker=marker1, s=150.0, facecolor=face,
                           edgecolor=edge, flipped=flipped, frame=frame,
                           **self._kwargs)
        point2 = sky_point(x, y, ax, marker=marker2, s=75.0, facecolor=face,
                           edgecolor=edge, flipped=flipped, frame=frame,
                           **self._kwargs)
        return [point1, point2]

    def __repr__(self):
        s = '<Sun: alpha={};\n'.format(self.alpha)
        s += '      size={}>'.format(self.size)
        return s


class DetectorPointing(SkyCircle):
    """Plot a detector pointing on the sky in equatorial, galactic, or 
    spacecraft coordinates.

    Parameters:
        x (float): The azimuthal coordinate, in degrees
        y (float): The polar coordinate, in degrees
        radius (float):  The radius of the circle, in degrees
        ax (:class:`matplotlib.axes`): The axis on which to plot
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', 'spacecraft'.
                               Default is 'equatorial'
        face_color (str, optional): The color of the circle fill
        face_alpha (float, optional): The alpha opacity of the circle fill
        edge_color (str, optional): The color of the circle edge
        edge_alpha (float, optional): The alpha opacity of the circle edge
        color (str, optional): The color of the circle. If set, overrides 
                               ``face_color`` and ``edge_color``.
        alpha (float, optional): The alpha of the circle. If set, overrides 
                               ``face_alpha`` and ``edge_alpha``.
        fontsize (float, optional): The size of the detector label
        font_alpha (float, optional): The alpha opacity of the detector label
        **kwargs: Other plotting options
    """
    def __init__(self, x, y, radius, det, ax, flipped=True, frame='equatorial',
                 color='dimgray', alpha=None, edge_alpha=0.5, face_alpha=0.25, 
                 fontsize=10, font_alpha=0.8, font_color=None, **kwargs):
        
        super().__init__(x, y, radius, ax, color=color, alpha=alpha, 
                         face_alpha=face_alpha, edge_alpha=edge_alpha, 
                         frame=frame)
        self._det = det
        self._fontsize = fontsize
        self._fontalpha = font_alpha
        self._font_color = font_color
        if font_color is None:
            self._font_color = self._color
        self._annotate(x, y, ax, flipped, frame)

    @property
    def font_alpha(self):
        """(float): The alpha opacity of the detector label"""
        return self._fontalpha
    @font_alpha.setter
    def font_alpha(self, alpha):
        self._artists[-1].set_alpha(alpha)
        self._fontalpha = alpha

    @property
    def font_color(self):
        """(float): The color of the detector label"""
        return self._font_color
    @font_color.setter
    def font_color(self, color):
        self._artists[-1].set_color(color)
        self._font_color = color

    @property
    def fontsize(self):
        """(float): The size of the detector label"""
        return self._fontsize
    @fontsize.setter
    def fontsize(self, size):
        self._artists[-1].set_fontsize(size)
        self._fontsize = size

    def _annotate(self, x, y, ax, flipped, frame):
        theta = np.deg2rad(y)
        phi = np.deg2rad(180.0 - x)
        if frame == 'spacecraft':
            phi -= np.pi
            if phi < -np.pi:
                phi += 2 * np.pi
        elif frame == 'galactic':
            phi -= np.pi
            if phi < -np.pi:
                phi += 2 * np.pi
        else:
            pass            

        if not flipped:
            phi *= -1.0

        txt = ax.text(phi, theta, self._det, fontsize=self._fontsize,
                      ha='center', va='center', color=self._font_color,
                      alpha=self._fontalpha, **self._kwargs)
        self._artists.append(txt)

    def __repr__(self):
        spaces = ' '*19
        s = "<DetectorPointing: '{}';\n".format(self._det)
        s += '{0}face_color={1};\n'.format(spaces, self.face_color) 
        s += "{0}face_alpha={1};\n".format(spaces, self.face_alpha)
        s += "{0}edge_color={1};\n".format(spaces, self.edge_color)
        s += "{0}edge_alpha={1};\n".format(spaces, self.edge_alpha)
        s += "{0}linestyle='{1}';\n".format(spaces, self.linestyle)
        s += '{0}linewidth={1};\n'.format(spaces, self.linewidth)
        s += '{0}fontsize={1};\n'.format(spaces, self.fontsize)
        s += '{0}font_color={1};\n'.format(spaces, self.font_color)
        s += '{0}font_alpha={1}>'.format(spaces, self.font_alpha)
        return s
    

class GalacticPlane(PlotElement):
    """Plot the Galactic Plane in equatorial, galactic, or spacecraft 
    coordinates.

    Parameters:
        ax (:class:`matplotlib.axes`): The axis on which to plot
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following equatorial convention
        frame (str or :class:`~gdt.core.coords.SpacecraftFrame`, optional): 
            If a string, then can either be 'equatorial' or 'galactic'.
            Otherwise, it is the spacecraft frame definition. Defaults is 
            'equatorial'.
        inner_color (str, optional): The color of the inner line element
        outer_color (str, optional): The color of the outer line element
        line_alpha (float, optional): The alpha opacity of the line elements of 
                                      the Galactic Plane
        center_alpha (float, optional): The alpha opacity of the Galactic center
        color (str, optional): The color of the Galactic plane. Overrides
                               ``inner_color`` and ``outer_color``.
        alpha (float, optional): The opacity of the Galactic plane. Overrides
                                 ``line_alpha`` and ``center_alpha``.
        **kwargs: Other plotting options
    """
    def __init__(self, ax, flipped=True, frame='equatorial', color=None,
                 inner_color='black', outer_color='dimgray', alpha=None,
                 line_alpha=0.5, center_alpha=0.75, **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        # color and alpha act as setting the global color values for the object
        if color is not None:
            inner_color = color
            outer_color = color
        if alpha is not None:
            line_alpha = alpha
            center_alpha = alpha
        self._line_alpha = line_alpha
        self._center_alpha = center_alpha
        self._inner_color = inner_color
        self._outer_color = outer_color

        artists = self._create(ax, flipped, frame)
        self._artists = self._sanitize_artists(artists)

    @property
    def alpha(self):
        """(float): The alpha opacity value"""
        return self._alpha
    @alpha.setter
    def alpha(self, alpha):
        self.line_alpha = alpha
        self.center_alpha = alpha
        self._alpha = alpha

    @property
    def center_alpha(self):
        """(float): The alpha opacity of the Galactic plane"""
        return self._center_alpha
    @center_alpha.setter
    def center_alpha(self, alpha):
        [artist.set_alpha(alpha) for artist in self._artists \
         if artist.__class__.__name__ == 'PathCollection']
        self._center_alpha = alpha

    @property
    def color(self):
        """(str): The color of the Galactic plane"""
        return self._color
    @color.setter
    def color(self, color):
        self.inner_color = color
        self.outer_color = color
        self._color = color

    @property
    def inner_color(self):
        """(str): The color of the inner line element"""
        return self._inner_color
    @inner_color.setter
    def inner_color(self, color):
        [artist.set_color(color) for artist in self._artists \
         if artist.__class__.__name__ == 'Line2D' and \
         artist.get_color() == self._inner_color]
        self._artists[-1].set_facecolor(color)
        self._inner_color = color

    @property
    def line_alpha(self):
        """(float): The alpha opacity of the line elements"""
        return self._line_alpha
    @line_alpha.setter
    def line_alpha(self, alpha):
        [artist.set_alpha(alpha) for artist in self._artists \
         if artist.__class__.__name__ == 'Line2D']
        self._line_alpha = alpha

    @property
    def outer_color(self):
        """(str): The color of the outer line element"""
        return self._outer_color
    @outer_color.setter
    def outer_color(self, color):
        [artist.set_color(color) for artist in self._artists \
         if artist.__class__.__name__ == 'Line2D' and \
         artist.get_color() == self._outer_color]
        self._artists[-2].set_facecolor(color)
        self._outer_color = color

    def _create(self, ax, flipped, frame):
        refs = galactic_plane(ax, flipped=flipped, frame=frame,
                              outer_color=self._outer_color,
                              inner_color=self._inner_color,
                              line_alpha=self._line_alpha,
                              center_alpha=self._center_alpha, **self._kwargs)
        return refs

    def __repr__(self):
        spaces = ' ' * 16
        s = '<GalacticPlane: outer_color={};\n'.format(self.outer_color) 
        s += "{0}inner_color={1};\n".format(spaces, self.inner_color)
        s += "{0}line_alpha={1};\n".format(spaces, self.line_alpha)
        s += "{0}center_alpha={1};\n".format(spaces, self.center_alpha)
        return s


class SkyHeatmap(PlotElement):
    """Plot a heatmap on the sky. By default, the heatmap opacity will scale
    from 0 (fully transparent) to 1 (fully opaque) corresponding to the colormap
    value, creating an alpha gradient.  This behavior can be adjust by setting 
    ``alpha_min`` and ``alpha_max``.
    
    Parameters:
        x (np.array):The azimuthal coordinate array of the heatmap grid
        y (np.array): The polar coordinate array of the heatmap grid
        heatmap (np.array): The heatmap array, of shape (``x.size``, ``y.size``)
        ax (:class:`matplotlib.axes`): The axis on which to plot
        color (str, optional): The colormap of the heatmap. Default is 'RdPu'
        alpha_min (float, optional): The alpha opacity for the minimum color 
                                     value. Default is 0.
        alpha_max (float, optional): The alpha opacity for the maximum color 
                                     value. Default is 1.
        norm (:class:`matplotlib.colors.Normalize` or similar, optional):
            The normalization used to scale the colormapping to the heatmap 
            values. This can be initialized by the defined matplotlib 
            normalizations or a custom normalization. 
        flipped (bool, optional): If True, the azimuthal axis is flipped, 
                                  following equatorial convention
        frame (str, optional): Either 'equatorial', 'galactic', 'spacecraft'.
                               Default is 'equatorial'
        **kwargs: Other plotting options
    """
    def __init__(self, x, y, heatmap, ax, flipped=True, frame='equatorial',
                 color=GdtCmap('RdPu'), alpha=None, norm=None, **kwargs):

        # set the normalization, and ensure that it is scaled to the data range
        if norm is None:
            norm = PowerNorm(gamma=0.3)
        self._norm = norm
        self._norm.vmin = heatmap.min()
        self._norm.vmax = heatmap.max()
        self._heatmap = heatmap

        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(x, y, heatmap, ax, flipped, frame)
        self._artists = self._sanitize_artists(artists)

        # set the colormap
        self.color = color
        self.color.set_callback(self._artists[0].changed)

    @property
    def color(self):
        """(:class:`GdtCmap`): The colormap"""
        return self._color
    @color.setter
    def color(self, color):
        if not isinstance(color, GdtCmap):
            raise TypeError('color must be of type GdtCmap')
        self._color = color
        self._artists[0].set_cmap(color.cmap)

    @property
    def norm(self):
        """(:class:matplotlib.colors.Normalize or similar): 
        The colormap normalization""" 
        return self._norm
    @norm.setter
    def norm(self, norm):
        norm.vmin = self._heatmap.min()
        norm.vmax = self._heatmap.max()
        self._artists[0].set_norm(norm)
        self._norm = norm

    def _create(self, x, y, heatmap, ax, flipped, frame):
        refs = sky_heatmap(x, y, heatmap, ax, cmap=self._color.cmap,
                           norm=self._norm, flipped=flipped, frame=frame, 
                           **self._kwargs)
        ax.grid(True)
        return [refs]

    def __repr__(self):
        spaces = ' ' * 13
        s = "<SkyHeatmap: color='{0}';\n{1}".format(self.color.name, spaces) 
        s += 'norm={0}>'.format(self.norm.__class__.__name__)
        return s


class ModelData(PlotElement):
    """Plot a set of data with errors derived from the model variances.
    
    Parameters:
        x (np.array): The x coordinates of the data plot points
        y (np.array): The y coordinates of the data plot points
        xerr (np.array): The uncertainty of the x points
        yerr (np.array): The uncertainty of the y points
        ax (:class:`matplotlib.axes`): The axis on which to plot
        ulmask (np.array(dtype=bool), optional):
             A boolean array of the same size as x and y indicating which 
             points are upper limits (True) and which are data points (False). 
        color (str, optional): The color of the model
        alpha (float, optional): The alpha of the fill. Default is 1
        **kwargs: Other plotting options
    """
    def __init__(self, x, y, xerr, yerr, ax, ulmask=None, color='C0',
                 alpha=1.0, **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(x, y, xerr, yerr, ulmask, ax)
        self._artists = self._sanitize_artists(artists)

    def _create(self, x, y, xerr, yerr, ulmask, ax):
        return ax.errorbar(x, y, xerr=xerr, yerr=yerr, uplims=ulmask,
                           fmt='none',
                           color=self._color, alpha=self._alpha,
                           **self._kwargs)

    def __repr__(self):
        spaces = ' ' * 12
        s = "<ModelData: color='{0}';\n{1}".format(self.color, spaces)
        s += "alpha={0}>".format(self.alpha) 
        return s


class ModelSamples(PlotElement):
    """Plot a series of spectral model samples

    Parameters:
        x (np.array): The x coordinates of the model
        y_list (list of np.array): A list of arrays, where each element is a 
                                  single spectral sample, having length equal 
                                  to x
        ax (:class:`matplotlib.axes`): The axis on which to plot
        color (str, optional): The color of the samples
        alpha (float, optional): The alpha of the fill. Default is 1
        label (str, optional): The label for the samples
        **kwargs: Other plotting options
    """
    def __init__(self, x, y_list, ax, color='C0', alpha=1.0, label=None,
                 **kwargs):
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(x, y_list, ax)
        self._artists = self._sanitize_artists(artists)
        if label is not None:
            self._artists[0].set_label(label)

    @property
    def linestyle(self):
        """(str): The linestyle of the histogram. Default is '-'"""
        return self._artists[0].get_linestyle()
    @linestyle.setter
    def linestyle(self, val):
        for artist in self._artists:
            artist.set_linestyle(val)
            
    @property
    def linewidth(self):
        """(int): The line width of the histogram. Default is 1.5"""
        return self._artists[0].get_linewidth()
    @linewidth.setter
    def linewidth(self, val):
        for artist in self._artists:
            artist.set_linewidth(val)

    def _create(self, x, y_list, ax):
        return [ax.plot(x, y, color=self._color, alpha=self._alpha, 
                        **self._kwargs) for y in y_list]

    def __repr__(self):
        spaces = ' ' * 15
        s = "<ModelSamples: color='{}';\n".format(self.color) 
        s += "{0}alpha={1};\n{0}linestyle='{2}';\n".format(spaces, self.alpha,
                                                        self.linestyle)
        s += '{0}linewidth={1}>'.format(spaces, self.linewidth)
        return s


class PlotElementCollection(DataCollection):
    """A collection class for plot elements.  This differs from DataCollection
    in that plot elements can be accessed as attributes by their name. 
    For example:
    
        >>> collection = PlotElementCollection()
        >>> collection.include(my_obj, name='yeehaw')
        >>> collection.yeehaw
        <my_obj>
    
    """
    def __getattr__(self, name):
        return self.get_item(name)
    