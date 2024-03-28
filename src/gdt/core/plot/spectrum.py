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
from ..data_primitives import ChannelBins
from .defaults import *
from .lib import *
from .plot import GdtPlot, Histo, HistoErrorbars, HistoFilled, \
    SpectrumBackground

__all__ = ['Spectrum']

class Spectrum(GdtPlot):
    """Class for plotting count spectra and count spectra paraphernalia. 
    
    This class can plot differential count spectra using energy information,
    or can plot a count spectrum based on raw energy channel.
    
    Parameters:
        data (:class:`~gdt.core.data_primitives.EnergyBins`, optional): 
            The count spectrum data to plot
        background (:class:`~gdt.background.primitives.BackgroundSpectrum`, optional): 
            The background spectrum to plot
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`
    """
    def __init__(self, data=None, background=None, canvas=None, **kwargs):
        super().__init__(canvas=canvas, **kwargs)

        self._spec = None
        self._errorbars = None
        self._bkgd = None
        self._selections = []

        # initialize the plot axes, labels, ticks, and scales
        self._ax.set_xlabel('Energy (keV)', fontsize=PLOTFONTSIZE)
        self._ax.set_ylabel('Rate (count/s-keV)', fontsize=PLOTFONTSIZE)
        self._ax.xaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.set_xscale('log')
        self._ax.set_yscale('log')

        # plot data and/or background if set on init
        if data is not None:
            self.set_data(data)
        if background is not None:
            self.set_background(background)

    @property
    def background(self):
        """(:class:`~gdt.plot.plot.SpectrumBackground`): The count spectrum 
        background plot element"""
        return self._bkgd

    @property
    def errorbars(self):
        """(:class:`~gdt.plot.plot.HistoErrorbars`): The error bars plot element
        """
        return self._errorbars

    @property
    def selections(self):
        """(list of :class:`~gdt.plot.plot.HistoFilled`): The count spectrum 
        selection plot element"""
        return self._selections

    @property
    def spectrum(self):
        """(:class:`~gdt.plot.plot.Histo`): The count spectrum plot element"""
        return self._spec

    # mark FIXME: store selections in a collection
    def add_selection(self, data):
        """Add a selection to the plot.  This adds a new selection to a list
        of existing selections.
        
        Args:
            data (:class:`~gdt.core.data_primitives.EnergyBins`): 
                The count spectrum data selections to plot
        """
        color, alpha, fill_alpha, kwargs = self._selection_settings()
        select = HistoFilled(data, self._ax, color=color, alpha=alpha,
                             fill_alpha=fill_alpha, **kwargs)
        self._selections.append(select)

    def set_background(self, background):
        """Set the background plotting data. If a background already exists,
        this triggers a replot of the background.
        
        Args:
            background (:class:`~gdt.background.primitives.BackgroundSpectrum`):
                The background spectrum to plot
        """
        color, alpha, band_alpha, kwargs = self._bkgd_settings()
        self._bkgd = SpectrumBackground(background, self._ax, color=color,
                                        alpha=alpha, band_alpha=band_alpha,
                                        zorder=1000, **kwargs)

    def set_data(self, data):
        """Set the count spectrum plotting data. If a count spectrum already 
        exists, this triggers a replot of the count spectrum.
        
        If an EnergyBins object is used, this will plot a differential energy
        spectrum (count/s/keV), and if a ChannelBins object is used (i.e. no
        energy information), the count spectrum per channel is plotted.
        
        Args:
            data (:class:`~gdt.core.data_primitives.EnergyBins` or 
                  :class:`~gdt.core.data_primitives.ChannelBins`):  The data
        """
        spec_color, spec_alpha, spec_kwargs = self._spec_settings()
        self._spec = Histo(data, self._ax, color=spec_color, alpha=spec_alpha,
                           **spec_kwargs)
        eb_color, eb_alpha, eb_kwargs = self._eb_settings()
        self._errorbars = HistoErrorbars(data, self._ax, color=eb_color,
                                         alpha=eb_alpha, **eb_kwargs)

        self._ax.set_xlim(data.range)
        mask = (data.rates > 0.0)
        
        if isinstance(data, ChannelBins):
            self._ax.set_xlabel('Channel Number', fontsize=PLOTFONTSIZE)
            self._ax.set_ylabel('Rate (count/s)', fontsize=PLOTFONTSIZE)
            self._ax.set_xscale('linear')
            self._ax.set_ylim(0.9 * data.rates[mask].min(),
                              1.1 * data.rates.max())
        else:
            self._ax.set_ylim(0.9 * data.rates_per_kev[mask].min(),
                              1.1 * data.rates_per_kev.max())
    
    def remove_background(self):
        """Remove the background from the plot.
        """
        self._bkgd.remove()
        self._bkgd = None

    def remove_data(self):
        """Remove the count spectrum from the plot.
        """
        self._spec.remove()
        self._spec = None

    def remove_errorbars(self):
        """Remove the count spectrum errorbars from the plot.
        """
        self._errorbars.remove()
        self._errorbars = None

    def remove_selections(self):
        """Remove the selections from the plot.
        """
        [selection.remove() for selection in self._selections]
        self._selections = []

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

    def _eb_settings(self):
        """The default settings for the errorbars. If a lightcurve already
        exists, use its errorbars settings instead.
        """
        if self._errorbars is None:
            eb_color = DATA_ERROR_COLOR
            eb_alpha = DATA_ERROR_ALPHA
            eb_kwargs = {}
        else:
            eb_color = self._errorbars.color
            eb_alpha = self._errorbars.alpha
            eb_kwargs = self._errorbars._kwargs
        return (eb_color, eb_alpha, eb_kwargs)

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

    def _spec_settings(self):
        """The default settings for the count spectrum. If a count spectrum 
        already exists, use its settings instead.
        """
        if self._spec is None:
            spec_color = DATA_COLOR
            spec_alpha = None
            spec_kwargs = {}
        else:
            spec_color = self._spec.color
            spec_alpha = self._spec.alpha
            spec_kwargs = self._spec._kwargs
        return (spec_color, spec_alpha, spec_kwargs)

