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
from gdt.core.coords import SpacecraftFrame, Quaternion
from gdt.core.healpix import HealPixLocalization
from gdt.core.plot.sky import SkyPlot, EquatorialPlot
from gdt.core.detector import Detectors

this_dir = os.path.dirname(__file__)


class TestDetectors(Detectors):
    __test__ = False
    det0 = ('Det0', 0,  45.0 * u.deg,  45.0 * u.deg)
    det1 = ('Det1', 1, 270.0 * u.deg, 135.0 * u.deg)


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

        cls.hpx.frame = SpacecraftFrame(
            quaternion=Quaternion.from_xyz_w(xyz=[0.0, 1.0, 0.0], w=1.0),
            obsgeoloc=r.CartesianRepresentation(-6320675.5, -1513143.1, 2313154.5, unit='m'),
            obstime=Time('2017-08-17 12:41:00.249', format='iso', scale='utc'),
            detectors=TestDetectors)

        cls.image_file = os.path.join(this_dir, "test.png")

    def tearDown(self):
        try:
            os.remove(self.image_file)
        except:
            pass

    def test_clevel_plot(self):
        plot = EquatorialPlot()
        plot.add_localization(self.hpx, clevels=[0.90, 0.50], gradient=False)
        plt.savefig(self.image_file) 

    def test_gradient_plot(self):
        plot = EquatorialPlot()
        plot.add_localization(self.hpx, clevels=[], gradient=True, detectors=[])
        plt.savefig(self.image_file)

    def test_frame_plot(self):
        plot = EquatorialPlot()
        plot.add_frame(self.hpx.frame)
        plt.savefig(self.image_file) 
