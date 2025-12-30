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
import unittest

from gdt.core.data_primitives import TimeEnergyBins,TimeChannelBins, Gti, Ebounds, EventList
from gdt.core.phaii import Phaii
from gdt.core.tte import PhotonList
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial
from gdt.core.background.unbinned import NaivePoisson

this_dir = os.path.dirname(__file__)

class TestBinnedFitter(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        counts = [[ 0,  0,  2,  1,  2,  0,  0,  0],
                  [ 3, 16, 10, 13, 14,  4,  3,  3],
                  [ 3, 23, 26, 13,  8,  8,  5,  5],
                  [ 4, 21, 19, 16, 13,  2,  3,  4],
                  [ 4, 20, 17, 11, 15,  2,  1,  5],
                  [ 6, 20, 19, 11, 11,  1,  4,  4]]
        tstart = [0.0000, 0.0039, 0.0640, 0.1280, 0.1920, 0.2560]          
        tstop = [0.0039, 0.0640, 0.1280, 0.1920, 0.2560, 0.320]
        exposure = [0.0038, 0.0598, 0.0638, 0.0638, 0.0638, 0.0638]
        emin = [4.323754, 11.464164, 26.22962, 49.60019, 101.016815,
                290.46063, 538.1436, 997.2431]
        emax = [11.464164, 26.22962, 49.60019, 101.016815, 290.46063,
                538.1436, 997.2431, 2000.]
         
        data = TimeEnergyBins(counts, tstart, tstop, exposure, emin, emax)
        gti = Gti.from_list([(0.0000, 0.320)])
        cls.phaii = Phaii.from_data(data, gti=gti, trigger_time=356223561.133346)
        
        cls.fitter = BackgroundFitter.from_phaii(cls.phaii, Polynomial)
        cls.fitter.fit(order=1)
    
    def test_dof(self):
        self.assertEqual(self.fitter.dof.size, 8)

    def test_livetime(self):
        self.assertEqual(self.fitter.livetime, 0.3188)

    def test_method(self):
        self.assertEqual(self.fitter.method, 'Polynomial')

    def test_parameters(self):
        self.assertDictEqual(self.fitter.parameters, {'order': 1})

    def test_statistic(self):
        self.assertEqual(self.fitter.statistic.size, 8)
    
    def test_statistic_name(self):
        self.assertEqual(self.fitter.statistic_name, 'chisq')

    def test_type(self):
        self.assertEqual(self.fitter.type, 'binned')

    def test_interpolate_bins(self):
        rates = self.fitter.interpolate_bins(self.phaii.data.tstart, 
                                             self.phaii.data.tstop)
        self.assertTupleEqual(rates.size, (6, 8))

    def test_write(self):

        filepath = os.path.join(this_dir, 'test.bak')

        rates = self.fitter.interpolate_bins(self.phaii.data.tstart,
                                             self.phaii.data.tstop)
        bak = rates.to_bak(time_range=(0, 1))
        bak.write(directory=this_dir, filename=os.path.basename(filepath), overwrite=True)

        self.assertTrue(os.path.exists(filepath))
        os.remove(filepath)

    def test_errors(self):
        with self.assertRaises(TypeError):
            BackgroundFitter.from_phaii(0.0, Polynomial)
        
        with self.assertRaises(NameError):
            BackgroundFitter.from_phaii(self.phaii, UndefinedClass)
        
        with self.assertRaises(NotImplementedError):
            BackgroundFitter.from_phaii(self.phaii, Phaii)
        
        with self.assertRaises(TypeError):
            BackgroundFitter.from_phaii(self.phaii, Polynomial, time_ranges=0.0)

class TestBinnedFitterUncalibrated(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        counts = [[ 0,  0,  2,  1,  2,  0,  0,  0],
                  [ 3, 16, 10, 13, 14,  4,  3,  3],
                  [ 3, 23, 26, 13,  8,  8,  5,  5],
                  [ 4, 21, 19, 16, 13,  2,  3,  4],
                  [ 4, 20, 17, 11, 15,  2,  1,  5],
                  [ 6, 20, 19, 11, 11,  1,  4,  4]]
        tstart = [0.0000, 0.0039, 0.0640, 0.1280, 0.1920, 0.2560]          
        tstop = [0.0039, 0.0640, 0.1280, 0.1920, 0.2560, 0.320]
        exposure = [0.0038, 0.0598, 0.0638, 0.0638, 0.0638, 0.0638]
        chan_nums=[0, 1, 2, 3, 4, 5, 6, 7]
         
        data = TimeChannelBins(counts, tstart, tstop, exposure, chan_nums)
        gti = Gti.from_list([(0.0000, 0.320)])
        cls.phaii = Phaii.from_data(data, gti=gti, trigger_time=356223561.133346)
        
        cls.fitter = BackgroundFitter.from_phaii(cls.phaii, Polynomial)
        cls.fitter.fit(order=1)
    
    def test_dof(self):
        self.assertEqual(self.fitter.dof.size, 8)

    def test_livetime(self):
        self.assertEqual(self.fitter.livetime, 0.3188)

    def test_method(self):
        self.assertEqual(self.fitter.method, 'Polynomial')

    def test_parameters(self):
        self.assertDictEqual(self.fitter.parameters, {'order': 1})

    def test_statistic(self):
        self.assertEqual(self.fitter.statistic.size, 8)
    
    def test_statistic_name(self):
        self.assertEqual(self.fitter.statistic_name, 'chisq')

    def test_type(self):
        self.assertEqual(self.fitter.type, 'binned')

    def test_interpolate_bins(self):
        rates = self.fitter.interpolate_bins(self.phaii.data.tstart, 
                                             self.phaii.data.tstop)
        self.assertTupleEqual(rates.size, (6, 8))

    def test_write(self):
        with self.assertRaises(NotImplementedError):
            rates=self.fitter.interpolate_bins(self.phaii.data.tstart,
                                             self.phaii.data.tstop)
            _ = rates.to_bak()

    def test_errors(self):
        with self.assertRaises(TypeError):
            BackgroundFitter.from_phaii(0.0, Polynomial)
        
        with self.assertRaises(NameError):
            BackgroundFitter.from_phaii(self.phaii, UndefinedClass)
        
        with self.assertRaises(NotImplementedError):
            BackgroundFitter.from_phaii(self.phaii, Phaii)
        
        with self.assertRaises(TypeError):
            BackgroundFitter.from_phaii(self.phaii, Polynomial, time_ranges=0.0)


class TestUnbinnedFitter(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        times = [0.706, 1.640, 3.185, 3.512, 5.540, 
                 7.790, 9.602, 9.726, 10.45, 10.61]
        chans = [4, 1, 0, 4, 5, 0, 4, 0, 2, 0]
        ebounds = Ebounds.from_bounds([10.0, 20.0, 40.0, 80.0, 160.0, 320.0], 
                                      [20.0, 40.0, 80.0, 160.0, 320.0, 640.0])
        data = EventList(times=times, channels=chans, ebounds=ebounds)
        gti = Gti.from_list([(0.0000, 10.70)])
        cls.tte = PhotonList.from_data(data, gti=gti, 
                                       trigger_time=356223561.133346,
                                       event_deadtime=0.001, 
                                       overflow_deadtime=0.1)
        
        cls.fitter = BackgroundFitter.from_tte(cls.tte, NaivePoisson)
        cls.fitter.fit(window_width=1.0)

    def test_dof(self):
        self.assertIsNone(self.fitter.dof)

    def test_livetime(self):
        self.assertEqual(self.fitter.livetime, 9.904)

    def test_method(self):
        self.assertEqual(self.fitter.method, 'NaivePoisson')

    def test_parameters(self):
        self.assertDictEqual(self.fitter.parameters, {'window_width': 1.0})

    def test_statistic(self):
        self.assertIsNone(self.fitter.statistic)
    
    def test_statistic_name(self):
        self.assertIsNone(self.fitter.statistic_name)

    def test_type(self):
        self.assertEqual(self.fitter.type, 'unbinned')

    def test_zeros(self):
        tstart = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        tstop = np.array([2.0, 3.0, 4.0, 4.0, 6.0])
        rates, uncert = self.fitter._method.interpolate(tstart, tstop)
        self.assertListEqual(list(rates[:,3]), [0.0, 0.0, 0.0, 0.0, 0.0])

    def test_interpolate_times(self):
    
        # simulated Poisson rate of 1 count/sec
        times = np.random.exponential(1.0, size=100).cumsum()
        # random channel numbers
        chans = np.random.randint(0, 6, size=100)
        # channel-to-energy mapping
        ebounds = Ebounds.from_bounds([10.0, 20.0, 40.0, 80.0, 160.0, 320.0],
                                      [20.0, 40.0, 80.0, 160.0, 320.0, 640.0])
        data = EventList(times=times, channels=chans, ebounds=ebounds)
        gti = Gti.from_list([(0.0000, 10.70)])
        tte = PhotonList.from_data(data, gti=gti, 
                                   trigger_time=356223561.133346,
                                   event_deadtime=0.001, 
                                   overflow_deadtime=0.1)
                                   
        fitter = BackgroundFitter.from_tte(tte, NaivePoisson)
        fitter.fit(window_width=10.0)
                                
        rates = fitter.interpolate_times(tte.data.times)
        self.assertTupleEqual(rates.size, (100, 6))
        self.assertListEqual(rates.tstart.tolist(), times.tolist())
        self.assertListEqual(rates.tstop.tolist(), times.tolist())
        self.assertListEqual(rates.emin.tolist(), ebounds.low_edges())
        self.assertListEqual(rates.emax.tolist(), ebounds.high_edges())
        self.assertListEqual(rates.exposure.tolist(), [0]*100)

    def test_errors(self):
        with self.assertRaises(TypeError):
            BackgroundFitter.from_tte(0.0, NaivePoisson)
        
        with self.assertRaises(NameError):
            BackgroundFitter.from_tte(self.tte, UndefinedClass)
        
        with self.assertRaises(NotImplementedError):
            BackgroundFitter.from_tte(self.tte, PhotonList)
        
        with self.assertRaises(TypeError):
            BackgroundFitter.from_tte(self.tte, Polynomial, time_ranges=0.0)
     
    
if __name__ == '__main__':
    unittest.main()
      
        
