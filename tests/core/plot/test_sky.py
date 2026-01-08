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

from gdt.core.healpix import HealPixLocalization
from gdt.core.plot.sky import SkyPlot, EquatorialPlot

this_dir = os.path.dirname(__file__)

class TestEquatorialPlot(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        nside = 8
        npix = hp.nside2npix(nside)

        # test map with single pixel smoothed by 0.2 radians
        arr = np.zeros(npix, dtype=float)
        arr[int(0.5 * npix)] = 1
        arr = hp.smoothing(arr, sigma=0.2)

        cls.hpx = HealPixLocalization.from_data(arr)

    def test_clevel_plot(self):
        plot = EquatorialPlot()
        plot.add_localization(self.hpx, clevels=[0.90, 0.50], gradient=False, detectors=[])

    def test_gradient_plot(self):
        plot = EquatorialPlot()
        plot.add_localization(self.hpx, clevels=[], gradient=True, detectors=[])
