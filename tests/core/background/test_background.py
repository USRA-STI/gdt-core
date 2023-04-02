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

from gdt.core.background.binned import Polynomial
from gdt.core.background.unbinned import NaivePoisson
from gdt.core.background.primitives import BackgroundRates, BackgroundSpectrum
from gdt.core.phaii import Phaii

class TestPolynomialBackground(TestCase):
    counts = np.array([[78.0, 58.0, 40.0, 26.0, 14.0, 6.0, 2.0, 0.0, 2.0, 6.0],
                       [6.0, 2.0, 0.0, 2.0, 6.0, 14.0, 26.0, 40.0, 58.0, 78.0]]).T
    exposure = np.array([2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0])
    edges = np.array([-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    tstart = edges[:-1]
    tstop = edges[1:]
    
    vals = np.array([[39.0, 29.0, 20.0, 13.0, 7.0, 3.0, 1.0, 0.0, 1.0, 3.0], 
                     [3.0, 1.0, 0.0, 1.0, 3.0, 7.0, 13.0, 20.0, 29.0, 39.0]]).T
    
    def test_fit(self): 
        bkgd = Polynomial(self.counts, self.tstart, self.tstop, self.exposure)
        bkgd.fit(order=2)
        self.assertEqual(bkgd.statistic_name, 'chisq')
        self.assertCountEqual(bkgd.dof, np.array([6, 6]))
        self.assertAlmostEqual(bkgd.statistic[0], 0.23, delta=0.01)
        self.assertAlmostEqual(bkgd.statistic[1], 0.23, delta=0.01)
    
    def test_interpolate(self):
        bkgd = Polynomial(self.counts, self.tstart, self.tstop, self.exposure)
        bkgd.fit(order=2)
        rates, rate_uncert = bkgd.interpolate(self.tstart, self.tstop)
        for i in range(10):
            self.assertAlmostEqual(rates[i,0], self.vals[i,0], delta=0.5)
            self.assertAlmostEqual(rates[i,1], self.vals[i,1], delta=0.5)
    

class TestNaivePoissonBackground(TestCase):
    times = [np.array([0., 1.14, 1.22, 1.28, 1.76, 3.45, 4.29, 4.78, 4.75, 5.42,
                      5.97, 7.40, 7.61, 7.98, 8.10, 8.16, 10.18, 10.13, 
                      13.22, 14.03])]*2
    
    def test_fit(self):
        bkgd = NaivePoisson(self.times)
        bkgd.fit(window_width=5.0, fast=True)
        bkgd.fit(window_width=5.0, fast=False)
    
    def test_interpolate(self):
        bkgd = NaivePoisson(self.times)
        bkgd.fit(window_width=5.0, fast=True)
        x = np.linspace(0.0, 14.0, 15)
        rates, uncert = bkgd.interpolate(x[:-1], x[1:])
        for i in range(14):
            self.assertAlmostEqual(rates[i,0], 1.0, delta=2.)
            self.assertAlmostEqual(rates[i,1], 1.0, delta=2.)

        bkgd.fit(window_width=5.0, fast=False)
        x = np.linspace(0.0, 14.0, 15)
        rates, uncert = bkgd.interpolate(x[:-1], x[1:])
        for i in range(14):
            self.assertAlmostEqual(rates[i,0], 1.0, delta=2.)
            self.assertAlmostEqual(rates[i,1], 1.0, delta=2.)

if __name__ == '__main__':
    unittest.main()
      
        
        