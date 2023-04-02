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

from gdt.core.background.primitives import BackgroundSpectrum
from gdt.core.data_primitives import ResponseMatrix, EnergyBins
from gdt.core.response import Rsp
from gdt.core.pha import Pha, Bak
from gdt.core.phaii import Phaii
from gdt.core.simulate.pha import PhaSimulator
from gdt.core.spectra.functions import Band, Comptonized

NUM_SIM = 10

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
    

class TestPhaSimulator(TestCase):
    
    def setUp(self):
        rsp = create_rsp()
        bkgd = create_bkgd()
        band_params = (0.01, 300.0, -1.0, -2.8)
        self.sim = PhaSimulator(rsp, Band(), band_params, 0.256, bkgd, 
                                'Gaussian')
    
    def test_simulate_background(self):
        bkgds = self.sim.simulate_background(NUM_SIM)
        self.assertEqual(len(bkgds), NUM_SIM)
        self.assertIsInstance(bkgds[0], BackgroundSpectrum)

    def test_simulate_source(self):
        srcs = self.sim.simulate_source(NUM_SIM)
        self.assertEqual(len(srcs), NUM_SIM)
        self.assertIsInstance(srcs[0], EnergyBins)

    def test_simulate_sum(self):
        specs = self.sim.simulate_sum(NUM_SIM)
        self.assertEqual(len(specs), NUM_SIM)
        self.assertIsInstance(specs[0], EnergyBins)

        nsim = 10000
        ref_rates = [128.55, 198.50, 97.25, 28.87]
        specs = self.sim.simulate_sum(nsim)
        rates = np.zeros((4, nsim))
        for i in range(nsim):
            rates[:,i] = specs[i].rates
        rates_mean = rates.mean(axis=1)
        
        for i in range(4):
            self.assertAlmostEqual(rates_mean[i], ref_rates[i], delta=2.0)

    def test_to_bak(self):
        baks = self.sim.to_bak(NUM_SIM)
        self.assertEqual(len(baks), NUM_SIM)
        self.assertIsInstance(baks[0], Bak)
        self.assertListEqual(baks[0].gti.low_edges(), [0])
        self.assertListEqual(baks[0].gti.high_edges(), [0.256])

        baks = self.sim.to_bak(1, tstart=10.0, tstop=10.3)
        self.assertListEqual(baks[0].gti.low_edges(), [10.0])
        self.assertListEqual(baks[0].gti.high_edges(), [10.3])

    def test_to_pha(self):
        phas = self.sim.to_pha(NUM_SIM)
        self.assertEqual(len(phas), NUM_SIM)
        self.assertIsInstance(phas[0], Pha)
        self.assertListEqual(phas[0].gti.low_edges(), [0])
        self.assertListEqual(phas[0].gti.high_edges(), [0.256])

        phas = self.sim.to_pha(1, tstart=10.0, tstop=10.3)
        self.assertListEqual(phas[0].gti.low_edges(), [10.0])
        self.assertListEqual(phas[0].gti.high_edges(), [10.3])

    def test_to_phaii(self):
        phaii = self.sim.to_phaii(NUM_SIM)
        self.assertIsInstance(phaii, Phaii)
        self.assertListEqual(phaii.gti.low_edges(), [0])
        self.assertListEqual(phaii.gti.high_edges(), [2.56])

        phaii = self.sim.to_phaii(NUM_SIM, bin_width=0.3)
        self.assertIsInstance(phaii, Phaii)
        self.assertListEqual(phaii.gti.low_edges(), [0])
        self.assertListEqual(phaii.gti.high_edges(), [3.00])

    def test_set_background(self):
        bkgd = create_bkgd()
        bkgd._rates = bkgd._rates/2.0
        self.sim.set_background(bkgd, 'Gaussian')
        
        NSIM = 10000
        bkgds = self.sim.simulate_background(NSIM)
        rates = np.zeros((4, NSIM))
        for i in range(NSIM):
            rates[:,i] = bkgds[i].rates
        rates_mean = rates.mean(axis=1)
        
        for i in range(4):
            self.assertAlmostEqual(rates_mean[i], bkgd.rates[i], delta=1.0)        

        self.sim.set_background(bkgd, 'Poisson')
        bkgds = self.sim.simulate_background(NSIM)
        rates = np.zeros((4, NSIM))
        for i in range(NSIM):
            rates[:,i] = bkgds[i].rates
        rates_mean = rates.mean(axis=1)
        
        for i in range(4):
            self.assertAlmostEqual(rates_mean[i], bkgd.rates[i], delta=1.0)        

    def test_set_rsp(self):
        rsp = create_rsp()
        rsp = Rsp.from_data(rsp.drm, start_time=-10.0, stop_time=10.0)
        self.sim.set_rsp(rsp)
        self.assertEqual(self.sim._src_gen._rsp.tstart, -10.0)
        self.assertEqual(self.sim._src_gen._rsp.tstop, 10.0)

    def test_set_source(self):
        comp_params = (0.01, 150.0, -0.5)
        self.sim.set_source(Comptonized(), comp_params, 1.024)
        
        ref_rates = [35.18, 66.79, 30.56, 0.02]
        NSIM = 10000
        srcs = self.sim.simulate_source(NSIM)
        rates = np.zeros((4, NSIM))
        for i in range(NSIM):
            rates[:,i] = srcs[i].rates
        rates_mean = rates.mean(axis=1)
        
        for i in range(4):
            self.assertAlmostEqual(rates_mean[i], ref_rates[i], delta=0.5)        

    def test_errors(self):
        with self.assertRaises(TypeError):
            self.sim.set_background(0, 'Poisson')

        with self.assertRaises(ValueError):
            self.sim.set_background(create_bkgd(), 'NIL')
        
        with self.assertRaises(TypeError):
            self.sim.set_rsp(0)
        
        with self.assertRaises(ValueError):
            self.sim.to_phaii(1, bin_width=0.1)
        

if __name__ == '__main__':
    unittest.main()
      