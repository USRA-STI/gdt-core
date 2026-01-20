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
import astropy.coordinates.representation as r
import astropy.units as u

from gdt.core.plot.lightcurve import Lightcurve
from gdt.core.data_primitives import TimeEnergyBins, Gti, TimeBins
from gdt.core.phaii import Phaii
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial

from . import MyMixin


class TestLightcurve(MyMixin, unittest.TestCase):

    @classmethod
    def setUpClass(cls):
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
        cls.selections = phaii.to_lightcurve(time_range=(0.2, 0.3))

    def test_lightcurve(self):
        l = Lightcurve(data=self.rates, background=self.back_rates)
        l.background.alpha = 0.28
        l.background.linestyle = ":"
        l.background.linewidth = 5.0

        plt.savefig(self.image_file)

        self.assertEqual(l.background.alpha, 0.28)
        self.assertEqual(l.background.linestyle, ":")
        self.assertEqual(l.background.linewidth, 5.0)
        self.assertEqual(str(l.background)[:21], "<LightcurveBackground")

    def test_selection(self):
        l = Lightcurve(data=self.rates, background=self.back_rates)
        l.add_selection(self.selections)
        self.assertEqual(len(l.selections), 1)
        plt.savefig(self.image_file)

    def test_remove(self):
        l = Lightcurve(data=self.rates, background=self.back_rates)
        l.add_selection(self.selections)

        l.remove_errorbars()
        l.remove_data()
        l.remove_background()
        l.remove_selections()

        self.assertEqual(l.errorbars, None)
        self.assertEqual(l.lightcurve, None)
        self.assertEqual(l.background, None)
        self.assertEqual(len(l.selections), 0)

        plt.savefig(self.image_file)

    def test_set(self):
        l = Lightcurve(data=self.rates, background=self.back_rates)
        l.add_selection(self.selections)

        # set all values again to test changes to existing lightcurve
        l.set_data(self.rates)
        l.set_background(self.back_rates)
        l.add_selection(self.selections)

        plt.savefig(self.image_file)
