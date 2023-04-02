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
from .plot import GdtPlot, Histo, HistoErrorbars, HistoFilled, \
    LightcurveBackground
from .lib import *
from .defaults import *

__all__ = ['Lightcurve']

class Lightcurve(GdtPlot):
    """Class for plotting lightcurves and lightcurve paraphernalia.
    
    Parameters:
        data (:class:`~gdt.core.data_primitives.TimeBins`, optional): 
            The lightcurve data to plot
        background (:class:`~gdt.background.primitives.BackgroundRates`, optional): 
            The background rates to plot
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`
    """
    def __init__(self, data=None, background=None, canvas=None, ax=None,
                 **kwargs):
        super().__init__(canvas=canvas, ax=ax, **kwargs)

        self._lc = None
        self._errorbars = None
        self._bkgd = None
        self._selections = []

        # initialize the plot axes, labels, ticks, and scales
        self._ax.set_xlabel('Time (s)', fontsize=PLOTFONTSIZE)
        self._ax.set_ylabel('Count Rate (count/s)', fontsize=PLOTFONTSIZE)
        self._ax.xaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.set_xscale('linear')
        self._ax.set_yscale('linear')

        # plot data and/or background if set on init
        if data is not None:
            self.set_data(data)
            self._ax.set_xlim(data.range)

            self._ax.set_ylim(0.9 * np.min(data.rates),
                              1.1 * np.max(data.rates))
        if background is not None:
            self.set_background(background)

    @property
    def background(self):
        """(:class:`~gdt.plot.plot.LightcurveBackground`): The background plot 
        element"""
        return self._bkgd

    @property
    def errorbars(self):
        """(:class:`~gdt.plot.plot.HistoErrorbars`): The error bars plot element
        """
        return self._errorbars

    @property
    def lightcurve(self):
        """(:class:`~gdt.plot.plot.Histo`): The lightcurve plot element"""
        return self._lc

    @property
    def selections(self):
        """(list of :class:`~gdt.plot.plot.HistoFilled`): The list of selection 
        plot elements"""
        return self._selections

    # mark FIXME: store selections in a collection
    def add_selection(self, data):
        """Add a selection to the plot.  This adds a new selection to a list
        of existing selections.
        
        Args:
            data (:class:`~gdt.core.data_primitives.TimeBins`): 
                The lightcurve data selection to plot
        """
        color, alpha, fill_alpha, kwargs = self._selection_settings()
        select = HistoFilled(data, self._ax, color=color, alpha=alpha,
                             fill_alpha=fill_alpha, **kwargs)
        self._selections.append(select)

    def remove_background(self):
        """Remove the background from the plot.
        """
        self._bkgd.remove()
        self._bkgd = None

    def remove_data(self):
        """Remove the lightcurve from the plot.
        """
        self._lc.remove()
        self._lc = None

    def remove_errorbars(self):
        """Remove the lightcurve error bars from the plot.
        """
        self._errorbars.remove()
        self._errorbars = None

    def remove_selections(self):
        """Remove the selections from the plot.
        """
        [selection.remove() for selection in self._selections]
        self._selections = []

    def set_background(self, background):
        """Set the background plotting data. If a background already exists,
        this triggers a replot of the background.
        
        Args:
            background (:class:`~gdt.background.primitives.BackgroundRates`): 
                The background model to plot
        """
        color, alpha, band_alpha, kwargs = self._bkgd_settings()
        self._bkgd = LightcurveBackground(background, self._ax, color=color,
                                          alpha=alpha, band_alpha=band_alpha,
                                          zorder=1000, **kwargs)

    def set_data(self, data):
        """Set the lightcurve plotting data. If a lightcurve already exists,
        this triggers a replot of the lightcurve.
        
        Args:
            data (:class:`~gdt.core.data_primitives.TimeBins`): 
                The lightcurve data to plot
        """
        lc_color, lc_alpha, lc_kwargs = self._lc_settings()
        self._lc = Histo(data, self._ax, color=lc_color, alpha=lc_alpha,
                         **lc_kwargs)
        eb_color, eb_alpha, eb_kwargs = self._eb_settings()
        self._errorbars = HistoErrorbars(data, self._ax, color=eb_color,
                                         alpha=eb_alpha, **eb_kwargs)

    def _lc_settings(self):
        """The default settings for the lightcurve. If a lightcurve already
        exists, use its settings instead.
        """
        if self._lc is None:
            lc_color = DATA_COLOR
            lc_alpha = None
            lc_kwargs = {}
        else:
            lc_color = self._lc.color
            lc_alpha = self._lc.alpha
            lc_kwargs = self._lc._kwargs
        return (lc_color, lc_alpha, lc_kwargs)

    def _eb_settings(self):
        """The default settings for the errorbars. If a lightcurve already
        exists, use its errorbars settings instead.
        """
        if self._errorbars is None:
            eb_color = DATA_ERROR_COLOR
            eb_alpha = None
            eb_kwargs = {}
        else:
            eb_color = self._errorbars.color
            eb_alpha = self._errorbars.alpha
            eb_kwargs = self._errorbars._kwargs
        return (eb_color, eb_alpha, eb_kwargs)

    def _bkgd_settings(self):
        """The default settings for the background. If a background already
        exists, use its settings instead.
        """
        if self._bkgd is None:
            color = BKGD_COLOR
            alpha = BKGD_ALPHA
            band_alpha = BKGD_ERROR_ALPHA
            kwargs = {'linewidth': BKGD_WIDTH}
        else:
            color = self._bkgd.color
            alpha = self._bkgd.alpha
            band_alpha = self._bkgd.band_alpha
            kwargs = self._bkgd._kwargs
        return color, alpha, band_alpha, kwargs

    def _selection_settings(self):
        """The default settings for a selection. If a selection already
        exists, use its settings instead.
        """
        if len(self._selections) == 0:
            color = DATA_SELECTED_COLOR
            alpha = 1.0
            fill_alpha = DATA_SELECTED_ALPHA
            kwargs = {}
        else:
            color = self._selections[0].color
            alpha = self._selections[0].alpha
            fill_alpha = self._selections[0].fill_alpha
            kwargs = self._selections[0]._kwargs
        return color, alpha, fill_alpha, kwargs
