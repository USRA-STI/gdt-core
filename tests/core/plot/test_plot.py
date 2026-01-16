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

from matplotlib.colors import Normalize
from gdt.core.plot.plot import *
from gdt.core.data_primitives import TimeEnergyBins, Gti, TimeBins
from gdt.core.phaii import Phaii
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial
from gdt.core.plot.sky import EquatorialPlot

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

    def test_cmap(self):
        cmap = GdtCmap("BuPu")

        def callback():
            print("callback function")
        cmap.set_callback(callback)

        cmap.alpha_scale = 'log'
        cmap.alpha_min = 0
        cmap.alpha_max = 0.5
        cmap.name = "Blues"

        self.assertEqual(cmap.alpha_scale, 'log')
        self.assertEqual(cmap.alpha_min, 1e-10)
        self.assertEqual(cmap.alpha_max, 0.5)
        self.assertEqual(cmap.name, "Blues")
        self.assertEqual(str(cmap)[:8], "<GdtCmap")

        with self.assertRaises(ValueError):
            cmap.alpha_min = -1
        with self.assertRaises(ValueError):
            cmap.alpha_min = 0.7
        with self.assertRaises(ValueError):
            cmap.alpha_max = -1
        with self.assertRaises(ValueError):
            cmap.alpha_min = 0.1
            cmap.alpha_max = 0.05
        with self.assertRaises(ValueError):
            cmap.alpha_scale = 'test'
        with self.assertRaises(ValueError):
            cmap = GdtCmap("BuPu", alpha_min=1000)
        with self.assertRaises(RuntimeError):
            cmap.set_callback([])

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

        self.assertEqual(h.linestyle, ":")
        self.assertEqual(h.linewidth, 2.0)
        self.assertEqual(h.alpha, 0.53)
        self.assertEqual(str(h)[:12], "<HistoFilled")

    def test_heatmap(self):
        fig = plt.figure()
        x = [0, 1, 2]
        y = [0, 1, 2]
        arr = np.array([[0, 1e-6, 1e-5], [1e-6, 1e-5, 1e-6], [1e-5, 1e-6, 0]])

        h = Heatmap(x, y, arr, plt.gca())
        plt.savefig(self.image_file)

        h.norm = Normalize()
        plt.savefig(self.image_file)

        with self.assertRaises(TypeError):
            h.color = 1.0

        self.assertEqual(str(h)[:8], "<Heatmap")

    def test_sky_circle(self):
        plot = EquatorialPlot()
        s = SkyCircle(270, 30, 5, plot.ax, alpha=0.22)
        plt.savefig(self.image_file)

        self.assertEqual(s.alpha, 0.22)

        s.fill = False
        s.hatch = "o"
        s.color = "g"
        s.alpha = 0.47
        s.linewidth = 5.0
        s.linestyle = ":"

        self.assertFalse(s.fill)
        self.assertEqual(s.hatch, "o")
        self.assertEqual(s.linewidth, 5.0)
        self.assertEqual(s.linestyle, ":")
        for c in [s.color, s.face_color, s.edge_color]:
            self.assertEqual(c, "g")
        for a in [s.alpha, s.face_alpha, s.edge_alpha]:
            self.assertEqual(a, 0.47)
        self.assertEqual(str(s)[:10], "<SkyCircle")

    def test_sky_annulus(self):
        plot = EquatorialPlot()
        s = SkyAnnulus(270, 30, 25, 1, plot.ax)
        plt.savefig(self.image_file)

        s.hatch = "o"
        s.color = "g"
        s.linewidth = 5.0
        s.linestyle = ":"

        self.assertEqual(s.hatch, "o")
        self.assertEqual(s.color, "g")
        self.assertEqual(s.linewidth, 5.0)
        self.assertEqual(s.linestyle, ":")
        self.assertEqual(str(s)[:11], "<SkyAnnulus")

    def test_sky_polygon(self):
        plot = EquatorialPlot()
        x = np.array([180, 270, 270, 180, 180])
        y = np.array([0, 0, 50, 50, 0])
        s = SkyPolygon(x, y, plot.ax)
        plt.savefig(self.image_file)

        s.fill = False
        s.hatch = "o"
        s.color = "g"
        s.alpha = 0.47
        s.linewidth = 5.0
        s.linestyle = ":"

        self.assertFalse(s.fill)
        self.assertEqual(s.hatch, "o")
        self.assertEqual(s.linewidth, 5.0)
        self.assertEqual(s.linestyle, ":")
        for c in [s.color, s.face_color, s.edge_color]:
            self.assertEqual(c, "g")
        for a in [s.alpha, s.face_alpha, s.edge_alpha]:
            self.assertEqual(a, 0.47)
        self.assertEqual(str(s)[:11], "<SkyPolygon")

    def test_sky_points(self):
        plot = EquatorialPlot()
        x = np.array([180, 270, 270, 180])
        y = np.array([0, 0, 50, 50])
        s = SkyPoints(x, y, plot.ax)
        s.sizes = 15
        plt.savefig(self.image_file)

        self.assertEqual(s.sizes[0], 15)
        self.assertEqual(s.num_points, 4)
        self.assertEqual(str(s)[:10], "<SkyPoints")

    def test_galactic_plane(self):
        plot = EquatorialPlot()
        g = GalacticPlane(plot.ax, color='m', alpha=0.66)
        plt.savefig(self.image_file)

        for c in [g.color, g.inner_color, g.outer_color]:
            self.assertEqual(c, "m")
        for a in [g.alpha, g.center_alpha, g.line_alpha]:
            self.assertEqual(a, 0.66)
        self.assertEqual(str(g)[:14], "<GalacticPlane")

        # test setters
        g.color = "g"
        self.assertEqual(g.color, "g")
        g.alpha = 0.33
        self.assertEqual(g.alpha, 0.33)
