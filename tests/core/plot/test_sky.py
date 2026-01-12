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
from gdt.core.healpix import HealPixLocalization, HealPixEffectiveArea
from gdt.core.plot.sky import SkyPlot, EquatorialPlot, GalacticPlot, SpacecraftPlot
from gdt.core.detector import Detectors
from gdt.core.plot.plot import PlotElementCollection

this_dir = os.path.dirname(__file__)


class MyDetectors(Detectors):
    det0 = ('Det0', 0,  45.0 * u.deg,  45.0 * u.deg)
    det1 = ('Det1', 1, 270.0 * u.deg, 135.0 * u.deg)


class MySkyPlot(SkyPlot):
    _x_start = 0
    _y_start = 0


class MyMixin:
    image_file = os.path.join(this_dir, "test.png")

    def tearDown(self):
        plt.close('all')
        try:
            os.remove(self.image_file)
        except:
            pass


class TestSkyPlot(unittest.TestCase, MyMixin):

    def test_fontsize(self):
        plot = MySkyPlot()
        self.assertEqual(plot.fontsize, 10)
        plot.fontsize = 12
        self.assertEqual(plot.fontsize, 12)

    def test_text_color(self):
        plot = MySkyPlot()
        self.assertEqual(plot.text_color, 'black')
        plot.text_color = 'red'
        self.assertEqual(plot.text_color, 'red')

    def test_properties(self):
        plot = MySkyPlot()
        self.assertEqual(plot.sun, None)
        self.assertEqual(plot.earth, None)
        self.assertEqual(plot.effective_area, None)
        self.assertEqual(plot.galactic_plane, None)
        self.assertEqual(plot.loc_posterior, None)
        self.assertEqual(len(plot.detectors), 0)
        self.assertIsInstance(plot.detectors, PlotElementCollection)
        self.assertEqual(len(plot.loc_contours), 0)
        self.assertIsInstance(plot.loc_contours, PlotElementCollection)

    def test_add_frame(self):
        w = 3 * [1]
        xyz = 3 * [[0.0, 1.0, 0.0]]
        eic = [
            (-6320675.5, -1513143.1, 2313154.5),
            (-6320575.5, -1513043.1, 2313154.5),
            (-6320475.5, -1512943.1, 2313154.5)
        ]
        time_str = ['2017-08-17 12:41:00.249', '2017-08-17 12:42:00.249', '2017-08-17 12:43:00.249']
        frames = SpacecraftFrame(
            quaternion=Quaternion.from_xyz_w(xyz=xyz, w=w),
            obsgeoloc=r.CartesianRepresentation(eic, unit='m'),
            obstime=Time(time_str, format='iso', scale='utc'),
            detectors=MyDetectors)

        plot = MySkyPlot()
        plot.add_frame(frames, trigtime=frames[1].obstime, detectors=[], sun=False, earth=False, galactic_plane=False)

        with self.assertRaises(ValueError):
            plot.add_frame(frames)
        with self.assertRaises(TypeError):
            plot.add_frame(frames, trigtime=0.0)

    def test_heatmap(self):
        az = np.linspace(0.0, 360, 10)
        el = np.linspace(-90, 90, 10)

        x, y = np.meshgrid(az, el)
        arr = np.zeros_like(x)
        for i in range(10):
            arr[i, :] = i + 1
        arr /= arr.sum()

        plot = MySkyPlot()
        heatmap = plot.plot_heatmap(arr, x, y)
        plt.savefig(self.image_file)

    def test_effective_area_plot(self):
        hpx = HealPixEffectiveArea.from_cosine(45, 45, 100, nside=8)

        plot = MySkyPlot()
        plot.add_effective_area(hpx)
        plt.savefig(self.image_file)


class TestEquatorialPlot(unittest.TestCase, MyMixin):

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
            detectors=MyDetectors)

    def test_add_localization(self):
        plot1 = EquatorialPlot()
        plot1.add_localization(self.hpx, gradient=True)
        plt.savefig(self.image_file)

        plot2 = EquatorialPlot()
        plot2.add_localization(self.hpx, gradient=False, detectors='det0')
        plt.savefig(self.image_file)

        # backup frame
        frame = self.hpx.frame

        # create dummy frame that will trigger exception pass for sun & earth plots
        class DummyFrame:
            detectors = frame.detectors

        self.hpx.frame = DummyFrame()

        plot3 = EquatorialPlot()
        plot3.add_localization(self.hpx, gradient=False, detectors=[])

        # restore frame
        self.hpx.frame = frame

    def test_frame_plot(self):
        plot1 = EquatorialPlot()
        plot1.add_frame(self.hpx.frame, detectors='det1')
        plt.savefig(self.image_file)

        plot2 = EquatorialPlot()
        plot2.add_frame(self.hpx.frame)

        self.assertEqual(len(plot1.detectors), 1)
        self.assertEqual(len(plot2.detectors), 2)

    def test_effective_area_plot(self):
        hpx = HealPixEffectiveArea.from_cosine(45, 45, 100, nside=8)

        plot1 = EquatorialPlot()
        plot1.add_effective_area(hpx, frame=self.hpx.frame, sun=True, earth=True, galactic_plane=True)
        plt.savefig(self.image_file)

        plot2 = EquatorialPlot()
        plot2.add_effective_area(hpx, frame=self.hpx.frame, detectors='all')
        plt.savefig(self.image_file)

        plot3 = EquatorialPlot()
        plot3.add_effective_area(hpx, frame=self.hpx.frame, detectors='det1')
        plt.savefig(self.image_file)

        with self.assertRaises(ValueError):
            plot3.add_effective_area(hpx)


class TestGalacticPlot(unittest.TestCase, MyMixin):

    @classmethod
    def setUpClass(cls):
        cls.frame = SpacecraftFrame(
            quaternion=Quaternion.from_xyz_w(xyz=[0.0, 0.0, 1.0], w=1.0),
            obsgeoloc=r.CartesianRepresentation(-6320675.5, -1513143.1, 2313154.5, unit='m'),
            obstime=Time('2017-08-17 12:41:00.249', format='iso', scale='utc'),
            detectors=MyDetectors)

    def test_effective_area_plot(self):
        hpx = HealPixEffectiveArea.from_cosine(45, 45, 100, nside=8)

        plot = GalacticPlot()
        plot.add_effective_area(hpx, frame=self.frame, sun=True, earth=True, galactic_plane=True)
        plt.savefig(self.image_file)

        with self.assertRaises(ValueError):
            plot.add_effective_area(hpx)
