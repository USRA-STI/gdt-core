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
from unittest import TestCase
import unittest

from gdt.core.background.primitives import BackgroundSpectrum
from gdt.core.data_primitives import ResponseMatrix, EnergyBins, EventList
from gdt.core.response import Rsp
from gdt.core.tte import PhotonList
from gdt.core.simulate.profiles import tophat, linear
from gdt.core.simulate.tte import TteBackgroundSimulator, TteSourceSimulator
from gdt.core.spectra.functions import Band, Comptonized


def create_rsp():
    # 8 photon bins x 4 energy channels
    matrix = [[25.2,  0.0,  0.0,  0.0],
              [51.8, 54.9,  0.0,  0.0],
              [2.59, 82.0, 44.8,  0.0],
              [3.10, 11.6, 77.0, 0.13],
              [1.26, 6.21, 29.3, 14.6],
              [0.45, 3.46, 13.8, 9.98],
              [0.52, 4.39, 13.3, 3.93],
              [0.79, 7.14, 16.1, 3.92]]
    emin = [5.00, 15.8, 50.0, 158., 500., 1581, 5000, 15811]
    emax = [15.8, 50.0, 158., 500., 1581, 5000, 15811, 50000]
    chanlo = [4.60, 27.3, 102., 538.]
    chanhi = [27.3, 102., 538., 2000]
    drm = ResponseMatrix(matrix, emin, emax, chanlo, chanhi)

    tstart = 524666421.47
    tstop = 524666521.47
    trigtime = 524666471.47
    rsp = Rsp.from_data(drm, start_time=tstart, stop_time=tstop,
                        trigger_time=trigtime)
    return rsp

def create_bkgd():
    rates = [37.5, 57.5, 17.5, 27.5]
    rate_uncert = [1.896, 2.889, 0.919, 1.66]
    emin = [4.60, 27.3, 102., 538.]
    emax = [27.3, 102., 538., 2000]
    exposure = 0.256
    back_spec = BackgroundSpectrum(rates, rate_uncert, emin, emax, exposure)
    return back_spec
    

class TestTteBackgroundSimulator(TestCase):
    
    def setUp(self):
        bkgd = create_bkgd()
        linear_params = (1.0, 0.5)
        self.sim = TteBackgroundSimulator(bkgd, 'Gaussian', linear, linear_params)

    def test_simulate(self):
        ev = self.sim.simulate(0.0, 10.0)
        self.assertIsInstance(ev, EventList)
        self.assertGreaterEqual(ev.time_range[0], 0.0)
        self.assertLessEqual(ev.time_range[1], 10.001)
        
        # integrating 1+0.5x linear background from (0, 10) seconds gives a
        # factor of 35.  background spectrum has an energy-integrated rate of
        # 140 counts/s.  140 * 35 = 4900 counts, therefore we should expect
        # ~4900 events with some variance.
        self.assertAlmostEqual(ev.size, 4900, delta=200)

    def test_to_tte(self):
        tte = self.sim.to_tte(0.0, 10.0)
        self.assertIsInstance(tte, PhotonList)
        self.assertGreaterEqual(tte.gti.low_edges()[0], 0.0)
        self.assertLessEqual(tte.gti.high_edges()[0], 10.001)

    def test_set_background(self):
        bkgd = create_bkgd()
        bkgd._rates = bkgd._rates/2.0

        self.sim.set_background(bkgd, 'Gaussian')        
        ev = self.sim.simulate(0.0, 10.0)
        self.assertAlmostEqual(ev.size, 2450, delta=200)
        

        self.sim.set_background(bkgd, 'Poisson')        
        ev = self.sim.simulate(0.0, 10.0)
        self.assertAlmostEqual(ev.size, 2450, delta=200)

    def test_errors(self):
        with self.assertRaises(ValueError):
            TteBackgroundSimulator(create_bkgd(), 'Gaussian', linear, 
                                   (1.0, 0.5), sample_period=0.0)

        with self.assertRaises(ValueError):
            TteBackgroundSimulator(create_bkgd(), 'Gaussian', linear, (1.0, 0.5),
                                   deadtime=-1.0)
        
        with self.assertRaises(TypeError):
            self.sim.set_background(0, 'Poisson')
        
        with self.assertRaises(ValueError):
            self.sim.set_background(create_bkgd(), 'NIL')


class TestTteSourceSimulator(TestCase):
    
    def setUp(self):
        rsp = create_rsp()
        band_params = (0.01, 300.0, -1.0, -2.8)
        tophat_params = (0.01, 1.0, 3.0)
        self.sim = TteSourceSimulator(rsp, Band(), band_params, tophat, 
                                      tophat_params)

    def test_simulate(self):
        ev = self.sim.simulate(0.0, 10.0)
        self.assertIsInstance(ev, EventList)
        self.assertGreaterEqual(ev.time_range[0], 1.0)
        self.assertLessEqual(ev.time_range[1], 3.001)
        
        # source spectrum has an energy-integrated rate of ~313 counts/s.  
        # therefore we should expect ~626 events over 2 s with some variance.
        self.assertAlmostEqual(ev.size, 626, delta=100)

    def test_to_tte(self):
        tte = self.sim.to_tte(0.0, 10.0)
        self.assertIsInstance(tte, PhotonList)
        self.assertGreaterEqual(tte.gti.low_edges()[0], 1.0)
        self.assertLessEqual(tte.gti.high_edges()[0], 3.001)

    def test_set_response(self):
        rsp = create_rsp()
        rsp = Rsp.from_data(rsp.drm, start_time=-10.0, stop_time=10.0)
        self.sim.set_response(rsp)
        self.assertEqual(self.sim._spec_gen._rsp.tstart, -10.0)
        self.assertEqual(self.sim._spec_gen._rsp.tstop, 10.0)

    def test_set_spectrum(self):
        comp_params = (0.01, 150.0, -0.5)
        self.sim.set_spectrum(Comptonized(), comp_params)
        ev = self.sim.simulate(0.0, 10.0)
        self.assertAlmostEqual(ev.size, 266, delta=50)
    
    def test_set_time_profile(self):
        self.sim.set_time_profile(tophat, (0.01, 1.0, 6.0))
        ev = self.sim.simulate(0.0, 10.0)
        self.assertAlmostEqual(ev.size, 1565, delta=200)

    def test_errors(self):
        with self.assertRaises(ValueError):
            TteSourceSimulator(create_rsp(), Band(), (0.01, 300.0, -1.0, -2.8), 
                               tophat, (0.01, 1.0, 3.0), sample_period=0.0)

        with self.assertRaises(ValueError):
            TteSourceSimulator(create_rsp(), Band(), (0.01, 300.0, -1.0, -2.8), 
                               tophat, (0.01, 1.0, 3.0), deadtime=-1.0)
        
        with self.assertRaises(TypeError):
            self.sim.set_response(0)


if __name__ == '__main__':
    unittest.main()
      