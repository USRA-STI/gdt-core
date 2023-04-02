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

from gdt.core.background.primitives import BackgroundRates, BackgroundSpectrum


class TestBackgroundRates(TestCase):

    def setUp(self):
        rates = [[30.0, 50.0, 10.0], [35.0, 55.0, 15.0], 
                 [40.0, 60.0, 20.0], [45.0, 65.0, 25.0]]
        rate_uncert = [[3.0, 5.0, 1.0], [3.5, 5.5, 1.5], 
                       [4.0, 6.0, 2.0], [4.5, 6.5, 2.5]]
        tstart = [0.0, 1.0, 3.0, 4.0]
        tstop = [1.0, 2.0, 4.0, 5.0]
        exposure = [2.0] * 4
        emin = [10.0, 50.0, 300.0]
        emax = [50.0, 150., 500.0]
        self.rates = BackgroundRates(rates, rate_uncert, tstart, tstop, emin,
                                     emax, exposure=exposure)
    
    def test_counts(self):
        counts = [[60, 100, 20], [70, 110, 30], [80, 120, 40], [90, 130, 50]]
        for i in range(4):
            self.assertListEqual(self.rates.counts[i,:].tolist(), counts[i])

    def test_count_uncertainty(self):
        uncert = [[6, 10, 2], [7, 11, 3], [8, 12, 4], [9, 13, 5]]
        for i in range(4):
            self.assertListEqual(self.rates.count_uncertainty[i,:].tolist(), 
                                 uncert[i])
    
    def test_rates(self):
        rates = [[30.0, 50.0, 10.0], [35.0, 55.0, 15.0], 
                 [40.0, 60.0, 20.0], [45.0, 65.0, 25.0]]
        for i in range(4):
            self.assertListEqual(self.rates.rates[i,:].tolist(), rates[i])
    
    def test_integrate_energy(self):
        # full range
        back_lc = self.rates.integrate_energy()
        self.assertListEqual(back_lc.rates.tolist(), [90.0, 105.0, 120.0, 135.0])
        uncert = [5.916, 6.690, 7.483, 8.292]
        for i in range(4):
            self.assertAlmostEqual(back_lc.rate_uncertainty[i], uncert[i], places=3)
        self.assertListEqual(back_lc.tstart.tolist(), self.rates.tstart.tolist())
        self.assertListEqual(back_lc.tstop.tolist(), self.rates.tstop.tolist())
        self.assertListEqual(back_lc.emin.tolist(), [10.0])
        self.assertListEqual(back_lc.emax.tolist(), [500.0])
        self.assertListEqual(back_lc.exposure.tolist(), self.rates.exposure.tolist())

        # set emin
        back_lc = self.rates.integrate_energy(emin=100.0)
        self.assertListEqual(back_lc.rates.tolist(), [60.0, 70.0, 80.0, 90.0])
        uncert = [5.099, 5.701, 6.325, 6.964]
        for i in range(4):
            self.assertAlmostEqual(back_lc.rate_uncertainty[i], uncert[i], places=3)
        self.assertListEqual(back_lc.emin.tolist(), [50.0])
        self.assertListEqual(back_lc.emax.tolist(), [500.0])

        # set emax
        back_lc = self.rates.integrate_energy(emax=100.0)
        self.assertListEqual(back_lc.rates.tolist(), [80.0, 90.0, 100.0, 110.0])
        uncert = [5.831, 6.519, 7.211, 7.906]
        for i in range(4):
            self.assertAlmostEqual(back_lc.rate_uncertainty[i], uncert[i], places=3)
        self.assertListEqual(back_lc.emin.tolist(), [10.0])
        self.assertListEqual(back_lc.emax.tolist(), [150.0])

        # set both emin and emax
        back_lc = self.rates.integrate_energy(emin=75., emax=125.)
        self.assertListEqual(back_lc.rates.tolist(), [50.0, 55.0, 60.0, 65.0])
        uncert = [5.0, 5.5, 6.0, 6.5]
        for i in range(4):
            self.assertAlmostEqual(back_lc.rate_uncertainty[i], uncert[i], places=3)
        self.assertListEqual(back_lc.emin.tolist(), [50.0])
        self.assertListEqual(back_lc.emax.tolist(), [150.0])

    def test_integrate_time(self):
        
        # full range
        back_spec = self.rates.integrate_time()
        self.assertListEqual(back_spec.rates.tolist(), [37.5, 57.5, 17.5])
        uncert = [1.896, 2.889, 0.919]
        for i in range(3):
            self.assertAlmostEqual(back_spec.rate_uncertainty[i], uncert[i], places=3)
        self.assertListEqual(back_spec.lo_edges.tolist(), [10.0, 50.0, 300.0])
        self.assertListEqual(back_spec.hi_edges.tolist(), [50.0, 150.0, 500.0])
        self.assertListEqual(back_spec.exposure.tolist(), [8.0]*3)

        # set tstart
        back_spec = self.rates.integrate_time(tstart=1.5)
        self.assertListEqual(back_spec.rates.tolist(), [40.0, 60.0, 20.0])
        uncert = [2.321, 3.472, 1.179]
        for i in range(3):
            self.assertAlmostEqual(back_spec.rate_uncertainty[i], uncert[i], places=3)
        self.assertListEqual(back_spec.exposure.tolist(), [6.0]*3)

        # set tstop
        back_spec = self.rates.integrate_time(tstop=3.5)
        self.assertListEqual(back_spec.rates.tolist(), [35.0, 55.0, 15.0])
        uncert = [2.034, 3.184, 0.898]
        for i in range(3):
            self.assertAlmostEqual(back_spec.rate_uncertainty[i], uncert[i], places=3)
        self.assertListEqual(back_spec.exposure.tolist(), [6.0]*3)

        # set both tstart amd tstop
        back_spec = self.rates.integrate_time(tstart=1.5, tstop=3.5)
        self.assertListEqual(back_spec.rates.tolist(), [37.5, 57.5, 17.5])
        uncert = [2.658, 4.070, 1.250]
        for i in range(3):
            self.assertAlmostEqual(back_spec.rate_uncertainty[i], uncert[i], places=3)
        self.assertListEqual(back_spec.exposure.tolist(), [4.0]*3)

    def test_rebin_energy(self):
        with self.assertRaises(NotImplementedError):
            self.rates.rebin_energy(lambda x:x)       

    def test_rebin_time(self):
        with self.assertRaises(NotImplementedError):
            self.rates.rebin_time(lambda x:x)       

    def test_slice_energy(self):
        # middle slice       
        bins2 = self.rates.slice_energy(20.0, 125.0)
        self.assertTupleEqual(bins2.energy_range, (10.0, 150.0))
        
        # slice below lower boundary
        bins2 = self.rates.slice_energy(1.0, 125.0)
        self.assertTupleEqual(bins2.energy_range, (10.0, 150.0))

        # slice above upper boundary
        bins2 = self.rates.slice_energy(125.0, 1000.0)
        self.assertTupleEqual(bins2.energy_range, (50.0, 500.0))
        
        # slice covering full range
        bins2 = self.rates.slice_energy(1.0, 1000.0)
        self.assertTupleEqual(bins2.energy_range, (10.0, 500.0))

        # slice one bin
        bins2 = self.rates.slice_energy(70.0, 70.0)
        self.assertTupleEqual(bins2.energy_range, (50.0, 150.0))

        # slice fully outside range
        bins2 = self.rates.slice_energy(1000.0, 2000.0)
        self.assertIsNone(bins2.energy_range)

    def test_slice_time(self):
        # middle slice       
        bins2 = self.rates.slice_time(1.5, 3.5)
        self.assertTupleEqual(bins2.time_range, (1.0, 4.0))
        
        # slice below lower boundary
        bins2 = self.rates.slice_time(-1.0, 3.5)
        self.assertTupleEqual(bins2.time_range, (0.0, 4.0))

        # slice above upper boundary
        bins2 = self.rates.slice_time(1.5, 10.0)
        self.assertTupleEqual(bins2.time_range, (1.0, 5.0))
        
        # slice covering full range
        bins2 = self.rates.slice_time(-1.0, 10.0)
        self.assertTupleEqual(bins2.time_range, (0.0, 5.0))

        # slice one bin
        bins2 = self.rates.slice_time(1.5, 1.5)
        self.assertTupleEqual(bins2.time_range, (1.0, 2.0))

        # slice fully outside range
        bins2 = self.rates.slice_time(10.0, 20.0)
        self.assertIsNone(bins2.time_range)

    def test_to_bak(self):
        bak = self.rates.to_bak()
        self.assertListEqual(bak.data.rates.tolist(), [37.5, 57.5, 17.5])
        uncert = [1.896, 2.889, 0.919]
        for i in range(3):
            self.assertAlmostEqual(bak.data.rate_uncertainty[i], uncert[i], places=3)
        self.assertListEqual(bak.data.lo_edges.tolist(), [10.0, 50.0, 300.0])
        self.assertListEqual(bak.data.hi_edges.tolist(), [50.0, 150.0, 500.0])
        self.assertListEqual(bak.data.exposure.tolist(), [8.0]*3)
        self.assertListEqual(bak.gti.low_edges(), [0.0])
        self.assertListEqual(bak.gti.high_edges(), [5.0])
    
    def test_merge_time(self):
        slice1 = self.rates.slice_time(0.5, 1.5)
        slice2 = self.rates.slice_time(3.5, 4.5)
        merged = BackgroundRates.merge_time([slice1, slice2])
        for i in range(4):
            self.assertListEqual(merged.rates[i,:].tolist(), 
                                 self.rates.rates[i,:].tolist())
        for i in range(4):
            self.assertListEqual(merged.rate_uncertainty[i,:].tolist(), 
                                 self.rates.rate_uncertainty[i,:].tolist())
        self.assertListEqual(merged.tstart.tolist(), self.rates.tstart.tolist())
        self.assertListEqual(merged.tstop.tolist(), self.rates.tstop.tolist())
        self.assertListEqual(merged.emin.tolist(), self.rates.emin.tolist())
        self.assertListEqual(merged.emax.tolist(), self.rates.emax.tolist())
        self.assertListEqual(merged.exposure.tolist(), self.rates.exposure.tolist())

    def test_sum_time(self):
        summed = BackgroundRates.sum_time([self.rates, self.rates])
        for i in range(4):
            self.assertListEqual(summed.rates[i,:].tolist(), 
                                 (2.0*self.rates.rates[i,:]).tolist())
        self.assertListEqual(summed.tstart.tolist(), self.rates.tstart.tolist())
        self.assertListEqual(summed.tstop.tolist(), self.rates.tstop.tolist())
        self.assertListEqual(summed.emin.tolist(), self.rates.emin.tolist())
        self.assertListEqual(summed.emax.tolist(), self.rates.emax.tolist())
        self.assertListEqual(summed.exposure.tolist(), self.rates.exposure.tolist())

        slice1 = self.rates.slice_time(1.5, 5.0)
        with self.assertRaises(AssertionError):
            BackgroundRates.sum_time([self.rates, slice1])
    
    def test_init_errors(self):
        with self.assertRaises(TypeError):
            BackgroundRates(0, self.rates.rate_uncertainty, self.rates.tstart,
                            self.rates.tstop, self.rates.emin, self.rates.emax,
                            exposure=self.rates.exposure)

        with self.assertRaises(TypeError):
            BackgroundRates(self.rates.rates.flatten(), 
                            self.rates.rate_uncertainty, self.rates.tstart,
                            self.rates.tstop, self.rates.emin, self.rates.emax,
                            exposure=self.rates.exposure)

        with self.assertRaises(TypeError):
            BackgroundRates(self.rates.rates, 0, self.rates.tstart,
                            self.rates.tstop, self.rates.emin, self.rates.emax,
                            exposure=self.rates.exposure)

        with self.assertRaises(TypeError):
            BackgroundRates(self.rates.rates, self.rates.rate_uncertainty.flatten(), 
                            self.rates.tstart, self.rates.tstop, 
                            self.rates.emin, self.rates.emax, 
                            exposure=self.rates.exposure)

        with self.assertRaises(TypeError):
            BackgroundRates(self.rates.rates, self.rates.rate_uncertainty.T, 
                            self.rates.tstart, self.rates.tstop, 
                            self.rates.emin, self.rates.emax, 
                            exposure=self.rates.exposure)


class TestBackgroundSpectrum(TestCase):
    
    def setUp(self):
        rates = [20.0, 50.0, 17.0, 3.0, 0.0]
        rate_uncertainty = [2.0, 5.0, 1.7, 0.3, 0.0]
        low = [10., 50., 180., 320., 500.]    
        high = [50., 180., 320., 500., 1000.]
        exposure = 2.0
        self.bins = BackgroundSpectrum(rates, rate_uncertainty, low, high, 
                                       exposure)
    
    def test_rates(self):
        self.assertListEqual(self.bins.rates.tolist(), 
                             [20.0, 50.0, 17.0, 3.0, 0.0])

    def test_rate_uncertainty(self):
        self.assertListEqual(self.bins.rate_uncertainty.tolist(), 
                             [2.0, 5.0, 1.7, 0.3, 0.0])

    def test_counts(self):
        self.assertListEqual(self.bins.counts.tolist(), 
                             [40.0, 100.0, 34.0, 6.0, 0.0])

    def test_count_uncertainty(self):
        self.assertListEqual(self.bins.count_uncertainty.tolist(), 
                             [4.0, 10.0, 3.4, 0.6, 0.0])

    def test_merge(self):
        with self.assertRaises(NotImplementedError):
            self.bins.merge([self.bins, self.bins])

    def test_rebin(self):
        with self.assertRaises(NotImplementedError):
            self.bins.rebin(lambda x: x)

    def test_slice(self):
        sliced = self.bins.slice(100.0, 350.0)
        self.assertListEqual(sliced.rates.tolist(), [50.0, 17.0, 3.0])     
        self.assertListEqual(sliced.rate_uncertainty.tolist(), [5.0, 1.7, 0.3])     
        self.assertListEqual(sliced.lo_edges.tolist(), [50., 180., 320.])     
        self.assertListEqual(sliced.hi_edges.tolist(), [180., 320., 500.])     
        self.assertListEqual(sliced.exposure.tolist(), [2.0]*3)     

    def test_sum(self):
        summed = BackgroundSpectrum.sum([self.bins, self.bins])
        self.assertListEqual(summed.rates.tolist(), [20.0, 50.0, 17.0, 3.0, 0.0])
        uncert = [1.414, 3.536, 1.202, 0.212, 0.0]
        for i in range(5):
            self.assertAlmostEqual(summed.rate_uncertainty[i], uncert[i], places=3)

        self.assertListEqual(summed.lo_edges.tolist(), self.bins.lo_edges.tolist())
        self.assertListEqual(summed.hi_edges.tolist(), self.bins.hi_edges.tolist())
        self.assertListEqual(summed.exposure.tolist(), [4.0]*5)

        slice1 = self.bins.slice(100.0, 350.0)
        with self.assertRaises(AssertionError):
            BackgroundSpectrum.sum([self.bins, slice1])

    def test_init_errors(self):
        with self.assertRaises(TypeError):
            BackgroundSpectrum(0.0, self.bins.rate_uncertainty, 
                               self.bins.lo_edges, self.bins.hi_edges, 
                               self.bins.exposure)

        with self.assertRaises(TypeError):
            BackgroundSpectrum(self.bins.rates, 0.0, 
                               self.bins.lo_edges, self.bins.hi_edges, 
                               self.bins.exposure)

        with self.assertRaises(TypeError):
            BackgroundSpectrum(self.bins.rates, self.bins.rate_uncertainty[1:], 
                               self.bins.lo_edges, self.bins.hi_edges, 
                               self.bins.exposure)

     
if __name__ == '__main__':
    unittest.main()
      
        
        