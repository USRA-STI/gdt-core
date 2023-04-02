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
import matplotlib.pyplot as  plt

from .defaults import *
from .plot import GdtPlot, Heatmap, EffectiveArea

__all__ = ['ResponsePlot', 'PhotonEffectiveArea', 'ChannelEffectiveArea']

class ResponsePlot(GdtPlot):
    """Class for plotting a response matrix.
    
    Parameters:
        drm (:class:`~gdt.core.data_primitives.ResponseMatrix`, optional): 
            The response object
        colorbar (bool, optional): If True, plot the colorbar for the 
                                   effective area scale. Default is True.
        ax (:class:`matplotlib.axes`, optional): 
            An existing axes object to plot to.  If not set, will create a new 
            axes object.
        multi (bool, optional):
            If True, plots a multiplot window showing the matrix and the 
            integrated effective area as a function of incident energy and 
            recorded energy.
        num_contours (int, optional): Number of contours to plot. Default is 100
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`
    """
    _background = 'black'

    def __init__(self, drm=None, colorbar=True, multi=False, canvas=None,
                 ax=None, num_contours=100, **kwargs):

        self._drm = None
        self._colorbar = colorbar
        self._multi = multi

        # do the multiplot
        if multi:
            self._colorbar = False
            ax, ax_x, ax_y = self._init_multiplot()
            self._p = PhotonEffectiveArea(drm=drm, canvas=canvas, ax=ax_x,
                                          **kwargs)
            self._c = ChannelEffectiveArea(drm=drm, canvas=canvas, ax=ax_y,
                                           orientation='horizontal', **kwargs)
            ax_x.get_xaxis().set_visible(False)
            ax_y.get_yaxis().set_visible(False)

        super().__init__(canvas=canvas, ax=ax, **kwargs)

        self._ax.set_facecolor(self._background)

        # initialize the plot axes, labels, ticks, and scales
        self._ax.set_xlabel('Photon Energy (keV)', fontsize=PLOTFONTSIZE)
        self._ax.set_ylabel('Channel Energy (keV)', fontsize=PLOTFONTSIZE)
        self._ax.xaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.set_xscale('log')
        self._ax.set_yscale('log')

        # plot data and/or background if set on init
        if drm is not None:
            self.set_response(drm, index=0, num_contours=num_contours)
    
    @property
    def drm(self):
        """(:class:`~gdt.plot.plot.Heatmap`): The matrix plot object"""
        return self._drm

    def set_response(self, drm, **kwargs):
        """Set the response data.
        
        Args:
            drm (:class:`~gdt.core.data_primitives.ResponseMatrix`): 
                The response object
            **kwargs: Arguments to pass to :class:`~gdt.plot.plot.Heatmap`
        """
        self._drm = Heatmap(drm.photon_bin_centroids, drm.channel_centroids,
                            drm.matrix.T, self.ax, colorbar=self._colorbar, 
                            **kwargs)

        self._ax.set_xlim(*drm.photon_bins.range)
        self._ax.set_ylim(*drm.ebounds.range)        

        # update the background color of the colorbar
        if self._colorbar:
            self._drm._artists[-1].patch.set_facecolor(self._background)
    
    def _init_multiplot(self):
        # initialize the multiplot

        # dimensions
        left, width = 0.12, 0.55
        bottom, height = 0.12, 0.55
        bottom_h = left_h = left + width
        matrix = [left, bottom, width, height]
        histx = [left, bottom_h, width, 0.40]
        histy = [left_h, bottom, 0.40, height]

        # create the plot axes
        main_ax = plt.axes(matrix)
        ax_x = plt.axes(histx)
        ax_y = plt.axes(histy)

        return (main_ax, ax_x, ax_y)


class PhotonEffectiveArea(GdtPlot):
    """Class for plotting the incident photon effective area.
    
    Parameters:
        drm (:class:`~gdt.core.data_primitives.ResponseMatrix`, optional): 
            The response object
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`
    """
    def __init__(self, drm=None, canvas=None, ax=None, **kwargs):
        super().__init__(canvas=canvas, ax=ax, **kwargs)
        self._drm = None

        # initialize the plot axes, labels, ticks, and scales
        self._ax.set_xlabel('Photon Energy (keV)', fontsize=PLOTFONTSIZE)
        self._ax.set_ylabel(r'Effective Area (cm$^2$)', fontsize=PLOTFONTSIZE)
        self._ax.xaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.set_xscale('log')
        self._ax.set_yscale('log')

        # plot data and/or background if set on init
        if drm is not None:
            self.set_response(drm)
    
    @property
    def drm(self):
        """(:class:`~gdt.plot.plot.EffectiveArea`): Effective area plot element
        """
        return self._drm
    
    def set_response(self, drm, **kwargs):
        """Set the response data.
        
        Args:
            drm (:class:`~gdt.core.data_primitives.ResponseMatrix`): 
                The response object
            **kwargs: Arguments to pass to :class:`~gdt.plot.plot.EffectiveArea`
        """
        _color, _alpha, _kwargs = self._settings()
        effarea = drm.photon_effective_area()
        self._drm = EffectiveArea(effarea, self._ax, color=_color,
                                   alpha=_alpha, **_kwargs)
        nz_mask = effarea.counts > 0
        self._ax.set_ylim(effarea.counts[nz_mask].min(), 
                          effarea.counts.max()*1.1)
        self._ax.set_xlim(drm.photon_bins[0].emin, drm.photon_bins[-1].emax)
        plt.draw()

    def _settings(self):
        """The default settings for the plot. If a plot already
        exists, use its settings instead.
        """
        if self.drm is None:
            _color = DATA_COLOR
            _alpha = None
            _kwargs = {}
        else:
            _color = self._data.color
            _alpha = self._data.alpha
            _kwargs = self._data._kwargs
        return (_color, _alpha, _kwargs)


class ChannelEffectiveArea(GdtPlot):
    """Class for plotting the recorded channel energy effective area.
    
    Parameters:
        drm (:class:`~gdt.core.data_primitives.ResponseMatrix`, optional): 
            The response object
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`
    """
    def __init__(self, drm=None, canvas=None, ax=None, orientation='vertical',
                 **kwargs):
        super().__init__(canvas=canvas, ax=ax, **kwargs)
        self._drm = None
        self._orientation = orientation

        # initialize the plot axes, labels, ticks, and scales
        if self._orientation == 'horizontal':
            self._ax.set_ylabel('Channel Energy (keV)', fontsize=PLOTFONTSIZE)
            self._ax.set_xlabel(r'Effective Area (cm$^2$)',
                                fontsize=PLOTFONTSIZE)
        else:
            self._ax.set_xlabel('Channel Energy (keV)', fontsize=PLOTFONTSIZE)
            self._ax.set_ylabel(r'Effective Area (cm$^2$)',
                                fontsize=PLOTFONTSIZE)
        self._ax.xaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.set_xscale('log')
        self._ax.set_yscale('log')

        # plot data and/or background if set on init
        if drm is not None:
            self.set_response(drm)

    @property
    def drm(self):
        """(:class:`~gdt.plot.plot.EffectiveArea`): Effective area plot element
        """
        return self._drm

    def set_response(self, drm, **kwargs):
        """Set the response data.
        
        Args:
            drm (:class:`~gdt.core.data_primitives.ResponseMatrix`): 
                The response object
            **kwargs: Arguments to pass to :class:`~gdt.plot.plot.EffectiveArea`
        """
        _color, _alpha, _kwargs = self._settings()
        effarea = drm.photon_effective_area()
        self._drm = EffectiveArea(effarea, self._ax, color=_color,
                                   alpha=_alpha, orientation=self._orientation,
                                   **_kwargs)

        xrange = (drm.ebounds[0].emin, drm.ebounds[-1].emax)
        nz_mask = effarea.counts > 0
        yrange = (effarea.counts[nz_mask].min(), effarea.counts.max()*1.1)
        if self._orientation == 'horizontal':
            self._ax.set_ylim(xrange)
            self._ax.set_xlim(yrange)
        else:
            self._ax.set_xlim(xrange)
            self._ax.set_ylim(yrange)

    def _settings(self):
        """The default settings for the plot. If a plor already
        exists, use its settings instead.
        """
        if self.drm is None:
            _color = DATA_COLOR
            _alpha = None
            _kwargs = {}
        else:
            _color = self._data.color
            _alpha = self._data.alpha
            _kwargs = self._data._kwargs
        return (_color, _alpha, _kwargs)
