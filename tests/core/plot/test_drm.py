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
import numpy as np
import unittest
import matplotlib.pyplot as plt

from gdt.core.plot.drm import ResponsePlot, PhotonEffectiveArea, ChannelEffectiveArea
from gdt.core.plot.plot import Heatmap
from gdt.core.data_primitives import ResponseMatrix

from . import MyMixin


class DrmMixin:

    @classmethod
    def setUpClass(cls):
        matrix = np.diag([0.1, 0.2, 0.3, 0.2, 0.1, 0.1])
        emin = [10., 20., 40., 80., 160., 320.]
        emax = [20., 40., 80., 160., 320., 640.]
        chanlo = [10.1, 20.1, 40.1, 80.1, 160.1, 320.1]
        chanhi = [20.1, 40.1, 80.1, 160.1, 320.1, 640.1]
        cls.drm = ResponseMatrix(matrix, emin, emax, chanlo, chanhi)


class TestDrmPlot(MyMixin, DrmMixin, unittest.TestCase):
    
    def test_plot(self):
        plot = ResponsePlot(drm=self.drm)
        plt.savefig(self.image_file)

        self.assertIsInstance(plot.drm, Heatmap)
        self.assertEqual(str(plot.drm)[:8], "<Heatmap")
        self.assertEqual(plot.drm.num_contours, 100)
        self.assertEqual(plot.drm.norm.vmin, 0.0)
        self.assertEqual(plot.drm.colorbar.vmax, 0.3)

    def test_multi(self):
        plot = ResponsePlot(drm=self.drm, multi=True)
        plt.savefig(self.image_file)

    def test_set(self):
        plot = ResponsePlot(drm=self.drm)
        plot.set_response(self.drm)
        plt.savefig(self.image_file)


class TestPhotonEffectiveAreaPlot(MyMixin, DrmMixin, unittest.TestCase):
    
    def test_plot(self):
        plot = PhotonEffectiveArea(drm=self.drm)
        plt.savefig(self.image_file)

        # set response again to check re-use of existing plot settings
        plot.set_response(self.drm)


class TestChannelEffectiveAreaPlot(MyMixin, DrmMixin, unittest.TestCase):
    
    def test_plot(self):
        plot = ChannelEffectiveArea(drm=self.drm)
        plot.drm.alpha = 0.28
        plot.drm.linestyle = ":"
        plot.drm.linewidth = 5.0

        plt.savefig(self.image_file)

        self.assertEqual(plot.drm.alpha, 0.28)
        self.assertEqual(plot.drm.linestyle, ":")
        self.assertEqual(plot.drm.linewidth, 5.0)
        self.assertEqual(str(plot.drm)[:14], "<EffectiveArea")

        # set response again to check re-use of existing plot settings
        plot.set_response(self.drm)
