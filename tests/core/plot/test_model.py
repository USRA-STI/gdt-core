#
#     Authors: William Cleveland (USRA),
#              Adam Goldstein (USRA) and
#              Daniel Kocevski (NASA)
#
#     Portions of the code are Copyright 2020 William Cleveland and
#     Adam Goldstein, Universities Space Research Association
#     All rights reserved.
#
#     Written for the Fermi Gamma-ray Burst Monitor (Fermi-GBM)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import os
import numpy as np
import healpy as hp
import unittest
import matplotlib.pyplot as plt
import astropy.coordinates.representation as r
import astropy.units as u

from astropy.time import Time
from gdt.core import data_path
from gdt.core.detector import Detectors
from gdt.core.coords import SpacecraftFrame, Quaternion
from gdt.core.spectra.fitting import SpectralFitterPgstat
from gdt.core.plot.model import ModelFit
from gdt.core.plot.plot import Histo, ModelSamples, PlotElementCollection

this_dir = os.path.dirname(__file__)


class MyMixin:
    image_file = os.path.join(this_dir, "test.png")

    @classmethod
    def setUpClass(cls):
        fit = data_path.joinpath('specfit.npz')
        cls.fitter = SpectralFitterPgstat.load(fit)

    def tearDown(self):
        plt.close('all')
        try:
            os.remove(self.image_file)
        except:
            pass


class TestModelFitPlot(MyMixin, unittest.TestCase):
    
    def test_plot(self):
        plot = ModelFit(fitter=self.fitter)
        self.assertEqual(plot.view, "counts")
        self.assertIsInstance(plot.count_data, PlotElementCollection)
        self.assertEqual(len(plot.count_data), 3)
        self.assertIsInstance(plot.count_models, PlotElementCollection)
        self.assertEqual(len(plot.count_models), 3)

        model_data = plot.count_data.get_item("n0")
        self.assertEqual(str(model_data)[:10], "<ModelData")

        plt.savefig(self.image_file)
        plt.ioff()

    def test_hide_residuals(self):
        plot = ModelFit(fitter=self.fitter)
        plot.hide_residuals()
        plt.savefig(self.image_file)
        plt.ioff()

    def test_show_residuals(self):
        plot = ModelFit(fitter=self.fitter, resid=False)
        plot.show_residuals(sigma=False)
        self.assertIsInstance(plot.residuals, PlotElementCollection)
        self.assertEqual(len(plot.residuals), 3)
        plt.savefig(self.image_file)
        plt.ioff()

    def test_photon_spectrum(self):
        plot = ModelFit(fitter=self.fitter, view="photon")
        self.assertIsInstance(plot.spectrum_model, PlotElementCollection)
        self.assertEqual(len(plot.spectrum_model), 3)

        model_samples = plot.spectrum_model.get_item("item1")
        model_samples.linestyle = ":"
        model_samples.linewidth = 5.0

        plt.savefig(self.image_file)

        self.assertEqual(model_samples.linestyle, ":")
        self.assertEqual(model_samples.linewidth, 5.0)
        self.assertEqual(str(model_samples)[:13], "<ModelSamples")

        plt.ioff()

    def test_energy_spectrum(self):
        plot = ModelFit(fitter=self.fitter)
        plot.energy_spectrum(plot_components=False, num_samples=10)
        plt.savefig(self.image_file)
        plt.ioff()

    def test_nufnu_spectrum(self):
        plot = ModelFit(fitter=self.fitter)
        plot.nufnu_spectrum(num_samples=10)
        plt.savefig(self.image_file)
        plt.ioff()

    def test_set(self):
        plot = ModelFit(fitter=self.fitter, view="test")
        plot._view = "photon"
        plot.set_fit(self.fitter)
        plot._view = "energy"
        plot.set_fit(self.fitter)
        plot._view = "nufnu"
        plot.set_fit(self.fitter)
        plt.ioff()
