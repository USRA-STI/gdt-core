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
from unittest import TestCase

from gdt.core.simulate.profiles import *


times = np.array([-10.0, 0.0, 10.0])


class TestTophat(TestCase):
    def test_tophat(self):
        params = (1.0, 0.0, 20.0)
        y = tophat(times, *params)
        self.assertCountEqual(y, np.array([0.0, 1.0, 1.0]))


class TestNorris(TestCase):
    def test_norris(self):
        params = (1.0, -1.0, 0.1, 2.0)
        y = norris(times, *params)
        true = np.array((0.0, 0.858, 0.006))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]


class TestConstant(TestCase):
    def test_constant(self):
        params = (1.0,)
        y = constant(times, *params)
        self.assertCountEqual(y, np.array([1.0, 1.0, 1.0]))


class TestLinear(TestCase):
    def test_linear(self):
        params = (1.0, -2.0,)
        y = linear(times, *params)
        self.assertCountEqual(y, np.array([21.0, 1.0, -19.0]))
    

class TestQuadratic(TestCase):
    def test_quadratic(self):
        params = (1.0, -2.0, 2.0)
        y = quadratic(times, *params)
        self.assertCountEqual(y, np.array([221.0, 1.0, 181.0]))
    

if __name__ == '__main__':
    unittest.main()
      
        
        