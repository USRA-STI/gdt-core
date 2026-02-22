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
import unittest
import matplotlib.pyplot as plt

from gdt.core.plot.lightcurve import Lightcurve
from gdt.core.data_primitives import TimeEnergyBins, Gti, ChannelBins
from gdt.core.phaii import Phaii
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial
from gdt.core.plot.spectrum import Spectrum

from . import ImageFileMixin


class TestLightcurve(ImageFileMixin, unittest.TestCase):

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
        cls.spec = phaii.to_spectrum()
        cls.back_spec = fitter.interpolate_bins(phaii.data.tstart, phaii.data.tstop).integrate_time(0, 0.320)
        cls.selections = phaii.to_lightcurve(time_range=(0.2, 0.3))

    def test_spectrum(self):
        s = Spectrum(data=self.spec, background=self.back_spec)
        s.background.alpha = 0.28
        s.background.linestyle = ":"
        s.background.linewidth = 5.0

        plt.savefig(self.image_file)

        self.assertEqual(s.background.alpha, 0.28)
        self.assertEqual(s.background.linestyle, ":")
        self.assertEqual(s.background.linewidth, 5.0)
        self.assertEqual(str(s.background)[:19], "<SpectrumBackground")

    def test_selection(self):
        s = Spectrum(data=self.spec, background=self.back_spec)
        s.add_selection(self.selections)
        self.assertEqual(len(s.selections), 1)
        plt.savefig(self.image_file)

    def test_remove(self):
        s = Spectrum(data=self.spec, background=self.back_spec)
        s.add_selection(self.selections)

        s.remove_errorbars()
        s.remove_data()
        s.remove_background()
        s.remove_selections()

        self.assertEqual(s.errorbars, None)
        self.assertEqual(s.spectrum, None)
        self.assertEqual(s.background, None)
        self.assertEqual(len(s.selections), 0)

        plt.savefig(self.image_file)

    def test_set(self):
        s = Spectrum(data=self.spec, background=self.back_spec)
        s.add_selection(self.selections)

        # set all values again to test changes to existing spectrum
        s.set_data(self.spec)
        s.set_background(self.back_spec)
        s.add_selection(self.selections)

        plt.savefig(self.image_file)

    def test_channel_bins(self):
        counts = [20, 50, 17, 3, 0, 3]
        chan_nums = [0, 1, 3, 4, 5, 6]
        exposure = 10.0
        channel_bins = ChannelBins.create(counts, chan_nums, exposure)

        s = Spectrum(data=channel_bins)

        plt.savefig(self.image_file)
