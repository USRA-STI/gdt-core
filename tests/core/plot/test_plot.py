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

from gdt.core.plot.plot import *
from gdt.core.data_primitives import TimeEnergyBins, Gti, TimeBins
from gdt.core.phaii import Phaii
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial

this_dir = os.path.dirname(__file__)


class TestPlot(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.image_file = os.path.join(this_dir, "test.png")

        counts = [98, 103, 100, 94, 105, 101, ]
        tstart = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        tstop = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        exposure = 6 * [1.0]

        cls.data = TimeBins(counts, tstart, tstop, exposure)

    def tearDown(self):
        plt.close('all')
        try:
            os.remove(self.image_file)
        except:
            pass

    def test_histo(self):
        fig = plt.figure()

        h = Histo(self.data, plt.gca(), label="test")
        h.linewidth = 2.0
        h.linestyle = ":"

        plt.savefig(self.image_file)

        self.assertEqual(h.linestyle, ":")
        self.assertEqual(h.linewidth, 2.0)
        self.assertEqual(str(h)[:6], "<Histo")

    def test_histo_errorbars(self):
        fig = plt.figure()

        h = HistoErrorbars(self.data, plt.gca(), label="test")
        h.linewidth = 2.0
        h.linestyle = ":"

        plt.savefig(self.image_file)

        self.assertEqual(h.linestyle, ":")
        self.assertEqual(h.linewidth, 2.0)
        self.assertEqual(str(h)[:15], "<HistoErrorbars")

    def test_histo_filled(self):
        fig = plt.figure()

        h = HistoFilled(self.data, plt.gca(), label="test")
        h.linewidth = 2.0
        h.linestyle = ":"
        h.alpha = 0.53

        plt.savefig(self.image_file)
        plt.show()

        self.assertEqual(h.linestyle, ":")
        self.assertEqual(h.linewidth, 2.0)
        self.assertEqual(h.alpha, 0.53)
        self.assertEqual(str(h)[:12], "<HistoFilled")
