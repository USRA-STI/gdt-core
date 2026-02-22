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

from astropy.time import Time
from gdt.core.coords import SpacecraftFrame, Quaternion
from gdt.core.plot.earthplot import EarthPlot
from gdt.core.geomagnetic import SouthAtlanticAnomaly
from gdt.core.plot.plot import EarthPoints, SAA, EarthLine

from . import ImageFileMixin, ExampleDetectors


class ExampleSaa(SouthAtlanticAnomaly):
    _latitude = [-30.0, 5.0, 5.0, -30.0, -30.0]
    _longitude = [-90.0, -90.0, 30.0, 30.0, -90.0]


class TestEarthPlot(ImageFileMixin, unittest.TestCase):

    def test_earth_plot(self):
        plot = EarthPlot()
        plt.savefig(self.image_file)
        self.assertNotEqual(plot.geoaxes, None)

    def test_frame(self):
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
            detectors=ExampleDetectors)

        plot1 = EarthPlot()
        plot1.add_spacecraft_frame(frames, trigtime=frames[1].obstime)
        plot1.standard_title()
        plt.savefig(self.image_file)

        plot2 = EarthPlot()
        plot2.add_spacecraft_frame(frames, trigtime=frames[1].obstime, icon=EarthPoints, marker='*')
        plot2.orbit.linestyle = ":"
        plot2.orbit.linewidth = 5.0
        plt.savefig(self.image_file)

        self.assertIsInstance(plot2.orbit, EarthLine)
        self.assertEqual(plot2.orbit.linestyle, ":")
        self.assertEqual(plot2.orbit.linewidth, 5.0)
        self.assertEqual(str(plot2.orbit)[:10], "<EarthLine")

        self.assertEqual(plot2.spacecraft.sizes[0], 10)
        self.assertEqual(str(plot2.spacecraft)[:12], "<EarthPoints")

    def test_saa(self):
        plot = EarthPlot(saa=ExampleSaa())
        plot.saa.linestyle = ":"
        plot.saa.linewidth = 5.0
        plot.saa.hatch = 'o'
        plot.saa.fill = False
        plt.savefig(self.image_file)

        arr = np.array([0])
        self.assertIsInstance(plot.saa, SAA)
        self.assertTrue(plot.saa.in_saa(arr, arr)[0])
        self.assertEqual(plot.saa.fill, False)
        self.assertEqual(plot.saa.hatch, "o")
        self.assertEqual(plot.saa.linestyle, ":")
        self.assertEqual(plot.saa.linewidth, 5.0)
        self.assertEqual(str(plot.saa)[:4], "<SAA")

    def test_errors(self):
        with self.assertRaises(TypeError):
            plot = EarthPlot(saa="test")
