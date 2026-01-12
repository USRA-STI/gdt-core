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

from gdt.core.plot.lightcurve import Lightcurve
from gdt.core.data_primitives import TimeEnergyBins, Gti, TimeBins
from gdt.core.phaii import Phaii
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial

this_dir = os.path.dirname(__file__)


class TestLightcurve(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.image_file = os.path.join(this_dir, "test.png")

        counts = [[ 0,  0,  2,  1,  2,  0,  0,  0],
                  [ 3, 16, 10, 13, 14,  4,  3,  3],
                  [ 3, 23, 26, 13,  8,  8,  5,  5],
                  [ 4, 21, 19, 16, 13,  2,  3,  4],
                  [ 4, 20, 17, 11, 15,  2,  1,  5],
                  [ 6, 20, 19, 11, 11,  1,  4,  4]]
        tstart = [0.0000, 0.0039, 0.0640, 0.1280, 0.1920, 0.2560]
        tstop = [0.0039, 0.0640, 0.1280, 0.1920, 0.2560, 0.320]
        exposure = [0.0038, 0.0598, 0.0638, 0.0638, 0.0638, 0.0638]
        emin = [4.323754, 11.464164, 26.22962, 49.60019, 101.016815,
                290.46063, 538.1436, 997.2431]
        emax = [11.464164, 26.22962, 49.60019, 101.016815, 290.46063,
                538.1436, 997.2431, 2000.]

        data = TimeEnergyBins(counts, tstart, tstop, exposure, emin, emax)
        gti = Gti.from_list([(0.0000, 0.320)])
        phaii = Phaii.from_data(data, gti=gti, trigger_time=356223561.133346)

        fitter = BackgroundFitter.from_phaii(phaii, Polynomial)
        fitter.fit(order=1)

        cls.phaii = phaii
        cls.rates = phaii.to_lightcurve()
        cls.back_rates = fitter.interpolate_bins(phaii.data.tstart, phaii.data.tstop)

    def tearDown(self):
        plt.close('all')
        try:
            os.remove(self.image_file)
        except:
            pass

    def test_properties(self):
        l = Lightcurve()
        self.assertEqual(l.errorbars, None)
        self.assertEqual(l.lightcurve, None)
        self.assertEqual(l.background, None)
        self.assertIsInstance(l.selections, list)
        self.assertEqual(len(l.selections), 0)

    def test_lightcurve(self):
        l = Lightcurve(data=self.rates, background=self.back_rates)
        plt.savefig(self.image_file)

    def test_selection(self):
        l = Lightcurve(data=self.rates, background=self.back_rates)
        l.add_selection(self.phaii.to_lightcurve(time_range=(0.2, 0.3)))
        self.assertEqual(len(l.selections), 1)
        plt.savefig(self.image_file)

    def test_remove(self):
        l = Lightcurve(data=self.rates, background=self.back_rates)
        l.add_selection(self.phaii.to_lightcurve(time_range=(0.2, 0.3)))

        l.remove_errorbars()
        l.remove_data()
        l.remove_background()
        l.remove_selections()

        self.assertEqual(l.errorbars, None)
        self.assertEqual(l.lightcurve, None)
        self.assertEqual(l.background, None)
        self.assertEqual(len(l.selections), 0)

        plt.savefig(self.image_file)
