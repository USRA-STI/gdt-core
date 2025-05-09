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
from gdt.core.data_primitives import ResponseMatrix
from gdt.core.response import Rsp
from gdt.core.simulate.generators import *
from gdt.core.spectra.functions import Band

NUM_SIM = 10000

rates = [37.5, 57.5, 17.5]
rate_uncert = [1.896, 2.889, 0.919]
emin = [10.0, 50.0, 150.0]
emax = [50.0, 150., 300.0]
exposure = 4.0
back_spec = BackgroundSpectrum(rates, rate_uncert, emin, emax, exposure)

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


class TestPoissonBackgroundGenerator(TestCase):
    
    def setUp(self):
        self.gen = PoissonBackgroundGenerator(back_spec)
    
    def test_validity(self):
        dev = next(self.gen)
        self.assertEqual(dev.size, 3)
        self.assertListEqual(dev.lo_edges.tolist(), emin)
        self.assertCountEqual(dev.hi_edges.tolist(), emax)
        self.assertCountEqual(dev.exposure.tolist(), [exposure]*3)
    
    def test_generate(self):
        rates_dev = np.zeros((3, NUM_SIM))
        for i in range(NUM_SIM):
            dev = next(self.gen)
            rates_dev[:,i] = dev.rates
        rates_mean = rates_dev.mean(axis=1)
        
        for i in range(3):
            self.assertAlmostEqual(rates_mean[i], rates[i], delta=0.5)        

    def test_set_rng(self):
        self.gen.set_rng(np.random.default_rng(seed=1))
        dev = next(self.gen)

        ref_rates = [37.5,  55.25, 19.75]
        for i in range(3):
            self.assertAlmostEqual(dev.rates[i], ref_rates[i])


class TestGaussianBackgroundGenerator(TestCase):
    
    def setUp(self):
        self.gen = GaussianBackgroundGenerator(back_spec)
    
    def test_validity(self):
        dev = next(self.gen)
        self.assertEqual(dev.size, 3)
        self.assertListEqual(dev.lo_edges.tolist(), emin)
        self.assertCountEqual(dev.hi_edges.tolist(), emax)
        self.assertCountEqual(dev.exposure.tolist(), [exposure]*3)
    
    def test_generate(self):
        rates_dev = np.zeros((3, NUM_SIM))
        for i in range(NUM_SIM):
            dev = next(self.gen)
            rates_dev[:,i] = dev.rates
        rates_mean = rates_dev.mean(axis=1)
        rates_std = rates_dev.std(axis=1)
        
        for i in range(3):
            self.assertAlmostEqual(rates_mean[i], rates[i], delta=0.5)        

        for i in range(3):
            self.assertAlmostEqual(rates_std[i], rate_uncert[i], delta=0.2)        
        

class TestSourceSpectrumGenerator(TestCase):
    
    def setUp(self):
        self.exposure = 0.256
        band_params = (0.01, 300.0, -1.0, -2.8)
        self.rsp = create_rsp()
        self.gen = SourceSpectrumGenerator(self.rsp, Band(), band_params, 
                                           self.exposure)
        self.rates = self.rsp.fold_spectrum(Band().fit_eval, band_params).rates
    
    def test_validity(self):
        dev = next(self.gen)
        self.assertEqual(dev.size, 4)
        self.assertListEqual(dev.lo_edges.tolist(), self.rsp.drm.ebounds.low_edges())
        self.assertCountEqual(dev.hi_edges.tolist(), self.rsp.drm.ebounds.high_edges())
        self.assertCountEqual(dev.exposure.tolist(), [self.exposure]*4)

    def test_generate(self):
        rates_dev = np.zeros((4, NUM_SIM))
        for i in range(NUM_SIM):
            dev = next(self.gen)
            rates_dev[:,i] = dev.rates
        rates_mean = rates_dev.mean(axis=1)
        
        for i in range(4):
            self.assertAlmostEqual(rates_mean[i], self.rates[i], delta=2.0)        


class TestVariablePoissonBackground(TestCase):

    def setUp(self):
        self.gen = VariablePoissonBackground(back_spec)
    
    def test_half(self):
        self.gen.amp = 0.5
        rates_dev = np.zeros((3, NUM_SIM))
        for i in range(NUM_SIM):
            dev = next(self.gen)
            rates_dev[:,i] = dev.rates
        rates_mean = rates_dev.mean(axis=1)
        
        for i in range(3):
            self.assertAlmostEqual(rates_mean[i], 0.5*rates[i], delta=0.5)

    def test_double(self):
        self.gen.amp = 2.0
        rates_dev = np.zeros((3, NUM_SIM))
        for i in range(NUM_SIM):
            dev = next(self.gen)
            rates_dev[:,i] = dev.rates
        rates_mean = rates_dev.mean(axis=1)
        
        for i in range(3):
            self.assertAlmostEqual(rates_mean[i], 2.0*rates[i], delta=0.5)
    
    def test_errors(self):
        with self.assertRaises(TypeError):
            self.gen.amp = 'foo'     

        with self.assertRaises(ValueError):
            self.gen.amp = -1.0 
        

class TestVariableGaussianBackground(TestCase):

    def setUp(self):
        self.gen = VariableGaussianBackground(back_spec)
    
    def test_half(self):
        self.gen.amp = 0.5
        rates_dev = np.zeros((3, NUM_SIM))
        for i in range(NUM_SIM):
            dev = next(self.gen)
            rates_dev[:,i] = dev.rates
        rates_mean = rates_dev.mean(axis=1)
        rates_std = rates_dev.std(axis=1)
        
        for i in range(3):
            self.assertAlmostEqual(rates_mean[i], 0.5*rates[i], delta=0.5)

        for i in range(3):
            self.assertAlmostEqual(rates_std[i], 0.5*rate_uncert[i], delta=0.2)        

    def test_double(self):
        self.gen.amp = 2.0
        rates_dev = np.zeros((3, NUM_SIM))
        for i in range(NUM_SIM):
            dev = next(self.gen)
            rates_dev[:,i] = dev.rates
        rates_mean = rates_dev.mean(axis=1)
        rates_std = rates_dev.std(axis=1)
        
        for i in range(3):
            self.assertAlmostEqual(rates_mean[i], 2.0*rates[i], delta=0.5)

        for i in range(3):
            self.assertAlmostEqual(rates_std[i], 2.0*rate_uncert[i], delta=0.2)        
    
    def test_errors(self):
        with self.assertRaises(TypeError):
            self.gen.amp = 'foo'     

        with self.assertRaises(ValueError):
            self.gen.amp = -1.0 


class TestSourceSpectrumGenerator(TestCase):
    
    def setUp(self):
        self.exposure = 0.256
        band_params = (0.01, 300.0, -1.0, -2.8)
        self.rsp = create_rsp()
        self.gen = VariableSourceSpectrumGenerator(self.rsp, Band(), band_params, 
                                           self.exposure)
        self.rates = self.rsp.fold_spectrum(Band().fit_eval, band_params).rates
    
    def test_half(self):
        self.gen.amp = 0.005
        rates_dev = np.zeros((4, NUM_SIM))
        for i in range(NUM_SIM):
            dev = next(self.gen)
            rates_dev[:,i] = dev.rates
        rates_mean = rates_dev.mean(axis=1)
        
        for i in range(4):
            self.assertAlmostEqual(rates_mean[i], 0.5*self.rates[i], delta=2.0)        

    def test_double(self):
        self.gen.amp = 0.02
        rates_dev = np.zeros((4, NUM_SIM))
        for i in range(NUM_SIM):
            dev = next(self.gen)
            rates_dev[:,i] = dev.rates
        rates_mean = rates_dev.mean(axis=1)
        
        for i in range(4):
            self.assertAlmostEqual(rates_mean[i], 2.0*self.rates[i], delta=2.0)        

    def test_errors(self):
        with self.assertRaises(TypeError):
            self.gen.amp = 'foo'     

        with self.assertRaises(ValueError):
            self.gen.amp = -1.0 
        

class TestEventSpectrumGenerator(TestCase):
    
    def setUp(self):
        self.count_spec = np.array([17, 19, 12,  1])
        self.gen = EventSpectrumGenerator(self.count_spec, 1e-6)
    
    def test_validity(self):
        times, chans = next(self.gen)
        self.assertEqual(times.size, self.count_spec.sum())
        self.assertEqual((chans == 0).sum(), 17)
        self.assertEqual((chans == 1).sum(), 19)
        self.assertEqual((chans == 2).sum(), 12)
        self.assertEqual((chans == 3).sum(), 1)
        self.assertGreaterEqual(times.min(), 0.0)
        self.assertLessEqual(times.max(), 1e-6)
        
    def test_update_spectrum(self):
        new_count_spec = np.array([46, 64, 42,  2])
        self.gen.spectrum = new_count_spec
        times, chans = next(self.gen)
        self.assertEqual(times.size, new_count_spec.sum())
        self.assertEqual((chans == 0).sum(), 46)
        self.assertEqual((chans == 1).sum(), 64)
        self.assertEqual((chans == 2).sum(), 42)
        self.assertEqual((chans == 3).sum(), 2)
        self.assertGreaterEqual(times.min(), 0.0)
        self.assertLessEqual(times.max(), 1e-6)
    
    def test_min_sep(self):
        gen = EventSpectrumGenerator(self.count_spec, 1e-4, min_sep=2.6e-6)
        times, chans = next(gen)
        sep = times[1:] - times[:-1]
        self.assertGreaterEqual(sep.min(), 2.6e-6)


if __name__ == '__main__':
    unittest.main()
      
        
