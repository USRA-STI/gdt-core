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
from .plot import GdtPlot, Histo, ModelData, ModelSamples
from .plot import PlotElementCollection as Collection
from .lib import *
from .defaults import *
import warnings
import matplotlib.pyplot as plt

__all__ = ['ModelFit']

class ModelFit(GdtPlot):
    """Class for plotting spectral fits.
    
    Parameters:
        fitter (:class:`~gdt.spectra.fitting.SpectralFitter`, optional): 
            The spectral fitter
        view (str, optional): The plot view, one of 'counts', 'photon', 
                              'energy' or 'nufnu'. Default is 'counts'
        resid (bool, optional): If True, plots the residuals in counts view. 
                                Default is True.
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`        
    """
    colors = '#7F3C8D,#11A579,#3969AC,#F2B701,#E73F74,#80BA5A,#E68310,#008695,#CF1C90,#f97b72,#4b4b8f,#A5AA99'.split(',')
    """(list): A list of default plotting colors to cycle through"""
    
    _min_y = 1e-10    

    def __init__(self, fitter=None, canvas=None, view='counts', resid=True,
                 interactive=True):
        
        warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
        
        self._figure, axes = plt.subplots(2, 1, sharex=True, sharey=False, 
                                          figsize=(5.7, 6.7), dpi=100, 
                                          gridspec_kw={'height_ratios': [3,1]})
        plt.subplots_adjust(hspace=0)
        self._ax = axes[0]
        self._resid_ax = axes[1]
        
        self._view = view
        self._fitter = None
        self._count_models = Collection()
        self._count_data = Collection()
        self._resids = Collection()
        self._spectrum_model = Collection()
                
                
        # plot data and/or background if set on init
        if fitter is not None:
            self.set_fit(fitter, resid=resid)
        if interactive:
            plt.ion()
    
    @property
    def count_data(self):
        """(:class:`~gdt.plot.plot.PlotElementCollection` of :class:`~gdt.plot.plot.ModelData`):
        The count data plot elements"""
        return self._count_data
    
    @property
    def count_models(self):
        """(:class:`~gdt.plot.plot.PlotElementCollection` of :class:`~gdt.plot.plot.Histo`): 
        The count model plot elements"""
        return self._count_models
    
    @property
    def spectrum_model(self):
        """(:class:`~gdt.plot.plot.PlotElementCollection` of :class:`~gdt.plot.plot.ModelSamples`):
        The model spectrum sample elements"""
        return self._spectrum_model
    
    @property
    def residuals(self):
        """(:class:`~gdt.plot.plot.PlotElementCollection` of :class:`~gdt.plot.plot.ModelData`):
        The fit residual plot elements"""
        return self._resids
    
    @property
    def view(self):
        """(str): The current plot view"""
        return self._view
        
    def count_spectrum(self):
        """Plot the count spectrum fit.
        """
        self._view = 'counts'
        self._ax.clear()
        
        model_counts = self._fitter.model_count_spectrum()
        energy, chanwidths, data_counts, data_counts_err, ulmasks = \
                                              self._fitter.data_count_spectrum()
        
        for i in range(self._fitter.num_sets):
            det = self._fitter.detectors[i]
            self._count_models.include(Histo(model_counts[i], self._ax, 
                                             edges_to_zero=False, 
                                             color=self.colors[i], alpha=1.0, 
                                             label=det), name=det)
            self._count_data.include(ModelData(energy[i], data_counts[i], 
                                               chanwidths[i], data_counts_err[i],
                                               self._ax, ulmask=ulmasks[i], 
                                               color=self.colors[i], 
                                               alpha=0.7, linewidth=0.9), 
                                               name=det)

        self._ax.set_ylabel(r'Rate [count s$^{-1}$ keV$^{-1}$]')
        self._set_view()
        self._ax.legend()
    
    def energy_spectrum(self, **kwargs):
        """Plot the energy spectrum model.
        
        Args:
            num_samples (int, optional): The number of sample spectra. 
                                         Default is 100.
            plot_components (bool, optional): Set to False to only plot the 
                                              overall model, not each component.
                                              Default is False.
        """
        self._view = 'energy'
        self._plot_spectral_model(**kwargs)
        self._ax.set_ylabel(r'Energy Flux [ph cm$^{-2}$ s$^{-1}$]', fontsize=PLOTFONTSIZE)

    def hide_residuals(self):
        """Hide the fit residuals.
        """
        try:
            self._figure.delaxes(self._resid_ax)
            self._ax.xaxis.set_tick_params(which='both', labelbottom=True)
            self._ax.set_xlabel('Energy (keV)', fontsize=PLOTFONTSIZE)
        except:
            print('Residuals already hidden')

    def nufnu_spectrum(self, **kwargs):
        """Plot the nuFnu spectrum model.
        
        Args:
            num_samples (int, optional): The number of sample spectra. 
                                         Default is 100.
            plot_components (bool, optional): Set to False to only plot the 
                                              overall model, not each component.
                                              Default is False.
        """
        self._view = 'nufnu'        
        self._plot_spectral_model(**kwargs)
        self._ax.set_ylabel(r'$\nu F_\nu$ [keV ph cm$^{-2}$ s$^{-1}$]', fontsize=PLOTFONTSIZE)
       
    def photon_spectrum(self, **kwargs):
        """Plot the photon spectrum model.
        
        Args:
            num_samples (int, optional): The number of sample spectra. 
                                         Default is 10.
            plot_components (bool, optional): Set to False to only plot the 
                                              overall model, not each component.
                                              Default is False.
        """
        self._view = 'photon'
        self._plot_spectral_model(**kwargs)
        self._ax.set_ylabel(r'Photon Flux [ph cm$^{-2}$ s$^{-1}$ keV$^{-1}$]', fontsize=PLOTFONTSIZE)

    def set_fit(self, fitter, resid=False):
        """Set the fitter. If a fitter already exists, this triggers a replot of 
        the fit.
        
        Args:
            fitter (:class:`~gdt.spectra.fitting.SpectralFitter`): 
                The spectral fitter for which a fit has been performed
            resid (bool, optional):  If True, plot the fit residuals
        """
        self._fitter = fitter
        
        if self._view == 'counts':
            self.count_spectrum()
            if resid:
                self.show_residuals()
            else:
                self.hide_residuals()
        
        elif self._view == 'photon':
            self.photon_spectrum()

        elif self._view == 'energy':
            self.energy_spectrum()

        elif self._view == 'nufnu':
            self.nufnu_spectrum()
        
        else:
            pass
           
    def show_residuals(self, sigma=True):
        """Show the fit residuals.
        
        Args:
            sigma (bool, optional): If True, plot the residuals in units of 
                                    model sigma, otherwise in units of counts. 
                                    Default is True.
        """
        # if we don't already have residuals axis
        if len(self._figure.axes) == 1:
            self._figure.add_axes(self._resid_ax)
        
        # get the residuals
        energy, chanwidths, resid, resid_err = self._fitter.residuals(sigma=sigma)
        
        # plot for each detector/dataset
        ymin, ymax = ([], [])
        for i in range(self._fitter.num_sets):
            det = self._fitter.detectors[i]
            self._resids.include(ModelData(energy[i], resid[i], chanwidths[i],
                                           resid_err[i], self._resid_ax, 
                                           color=self.colors[i], alpha=0.7, 
                                           linewidth=0.9), name=det)
            ymin.append((resid[i]-resid_err[i]).min())
            ymax.append((resid[i]+resid_err[i]).max())
        # the zero line
        self._resid_ax.axhline(0.0, color='black')
        self._resid_ax.set_xlabel('Energy [kev]', fontsize=PLOTFONTSIZE)
        
        if sigma:
            self._resid_ax.set_ylabel('Residuals [sigma]', fontsize=PLOTFONTSIZE)
        else:
            self._resid_ax.set_ylabel('Residuals [counts]', fontsize=PLOTFONTSIZE)
        
        # we have to set the y-axis range manually, because the y-axis
        # autoscale is broken (known issue) in matplotlib for this situation
        ymin = np.min(ymin)
        ymax = np.max(ymax)
        self._resid_ax.set_ylim((1.0-np.sign(ymin)*0.1)*ymin, 
                                (1.0+np.sign(ymax)*0.1)*ymax)
    
    def _plot_spectral_model(self, num_samples=100, plot_components=True):
        """Plot the spectral model by sampling from the Gaussian approximation
        to the parameters' posterior.
        
        Args:
            num_samples (int, optional): The number of sample spectra. 
                                         Default is 100.
        """
        # clean plot and hide residuals if any
        warnings.filterwarnings("ignore", category=UserWarning)
        self._ax.clear()
        self.hide_residuals()
        
        
        num_comp = self._fitter.num_components
        comps = self._fitter.function_components
        name = self._fitter.function_name
        
        self._spectrum_model = Collection()    
        
        # if the number of model components is > 1, plot each one
        if (num_comp > 1) and (plot_components):
            energies, samples = self._fitter.sample_spectrum(which=self._view, 
                                                             num_samples=num_samples, 
                                                             components=True)
            for i in range(num_comp):
                model = ModelSamples(energies, samples[:,i,:], self._ax,
                                     label=comps[i], color=self.colors[i+1],
                                     alpha=0.1, lw=0.3)
                self._spectrum_model.include(model)
                                            
            samples = samples.sum(axis=1)
        else:
            # or just plot the function
            energies, samples = self._fitter.sample_spectrum(which=self._view, 
                                                             num_samples=num_samples)
        y_max = samples.max(axis=(1,0))
        self._spectrum_model.include(ModelSamples(energies, samples, self._ax, 
                                        label=name, color=self.colors[0], 
                                        alpha=0.1, lw=0.3))
        self._set_view()
                
        # fix the alphas for the legend
        legend = self._ax.legend()
        for lh in legend.legendHandles:
            lh.set_alpha(1)
            lh.set_linewidth(1.0)

        if self._ax.get_ylim()[0] < self._min_y:
            self._ax.set_ylim(self._min_y, 10.0*y_max)

    def _set_view(self):
        """Set the view properties
        """
        self._ax.set_xlim(self._fitter.energy_range)
        self._ax.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self._ax.set_xscale('log')
        self._ax.set_yscale('log')
        self._ax.set_xlabel('Energy [kev]', fontsize=PLOTFONTSIZE)

