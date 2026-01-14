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
from gdt.core.detector import Detectors
from gdt.core.coords import SpacecraftFrame, Quaternion
from gdt.core.response import Rsp, Rsp2
from gdt.core.plot.drm import ResponsePlot, PhotonEffectiveArea, ChannelEffectiveArea
from gdt.core.plot.plot import Heatmap
from gdt.core.data_primitives import ResponseMatrix, Ebounds, EnergyBins, TimeBins

this_dir = os.path.dirname(__file__)


class MyMixin:
    image_file = os.path.join(this_dir, "test.png")

    @classmethod
    def setUpClass(cls):
        matrix = np.diag([0.1, 0.2, 0.3, 0.2, 0.1, 0.1])
        emin = [10., 20., 40., 80., 160., 320.]
        emax = [20., 40., 80., 160., 320., 640.]
        chanlo = [10.1, 20.1, 40.1, 80.1, 160.1, 320.1]
        chanhi = [20.1, 40.1, 80.1, 160.1, 320.1, 640.1]
        cls.drm = ResponseMatrix(matrix, emin, emax, chanlo, chanhi)

    def tearDown(self):
        plt.close('all')
        try:
            os.remove(self.image_file)
        except:
            pass


class TestDrmPlot(MyMixin, unittest.TestCase):
    
    def test_plot(self):
        plot = ResponsePlot(drm=self.drm)
        plt.savefig(self.image_file)
        self.assertIsInstance(plot.drm, Heatmap)

    def test_multi(self):
        plot = ResponsePlot(drm=self.drm, multi=True)
        plt.savefig(self.image_file)

    def test_set(self):
        plot = ResponsePlot(drm=self.drm)
        plot.set_response(self.drm)
        plt.savefig(self.image_file)


class TestPhotonEffectiveAreaPlot(MyMixin, unittest.TestCase):
    
    def test_plot(self):
        plot = PhotonEffectiveArea(drm=self.drm)
        plt.savefig(self.image_file)

        # set response again to check re-use of existing plot settings
        plot.set_response(self.drm)


class TestChannelEffectiveAreaPlot(MyMixin, unittest.TestCase):
    
    def test_plot(self):
        plot = ChannelEffectiveArea(drm=self.drm)
        plt.savefig(self.image_file)

        # set response again to check re-use of existing plot settings
        plot.set_response(self.drm)
