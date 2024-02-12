# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
# Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland and Adam Goldstein
#               Universities Space Research Association
#               Science and Technology Institute
#               https://sti.usra.edu
#
# Developed by: Daniel Kocevski
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#
import os
import unittest
from gdt.core.data_primitives import TimeEnergyBins, Gti, Ebounds, TimeBins, \
                                     TimeChannelBins, ChannelBins, EnergyBins
from gdt.core.headers import FileHeaders
from gdt.core.phaii import Phaii
from gdt.core.pha import Pha
from gdt.core.binning.binned import combine_by_factor

this_dir = os.path.dirname(__file__)

class TestPhaii(unittest.TestCase):
    
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

    def test_data(self):
        self.assertIsInstance(self.phaii.data, TimeEnergyBins)
    
    def test_ebounds(self):
        self.assertIsInstance(self.phaii.ebounds, Ebounds)
    
    def test_energy_range(self):
        self.assertTupleEqual(self.phaii.energy_range, (4.323754, 2000.))

    def test_filename(self):
        self.assertIsNone(self.phaii.filename)
    
    def test_gti(self):
        self.assertIsInstance(self.phaii.gti, Gti)
    
    def test_headers(self):
        self.assertIsNone(self.phaii.headers)

    def test_num_chans(self):
        self.assertEqual(self.phaii.num_chans, 8)
    
    def test_time_range(self):
        self.assertTupleEqual(self.phaii.time_range, (0.0, 0.320))
    
    def test_trigtime(self):
        self.assertEqual(self.phaii.trigtime, 356223561.133346)

    def test_get_exposure(self):
        # total exposure
        self.assertAlmostEqual(self.phaii.get_exposure(), 0.3188)
        
        # exposure of a time slice
        self.assertAlmostEqual(self.phaii.get_exposure(time_ranges=(0.0, 0.1)),
                               0.1274)

        # exposure of multiple time slices
        self.assertAlmostEqual(self.phaii.get_exposure(time_ranges=[(0.0, 0.1),
                                                                    (0.2, 0.3)]),
                               0.255)
    
    def test_rebin_energy(self):
        # full range
        rebinned_phaii = self.phaii.rebin_energy(combine_by_factor, 2)
        self.assertEqual(rebinned_phaii.num_chans, 4)
        
        # part of the range
        rebinned_phaii = self.phaii.rebin_energy(combine_by_factor, 2, 
                                                 energy_range=(25.0, 750.0))
        self.assertEqual(rebinned_phaii.num_chans, 5)

    def test_rebin_time(self):
        # full range
        rebinned_phaii = self.phaii.rebin_time(combine_by_factor, 2)
        self.assertEqual(rebinned_phaii.data.num_times, 3)

        # part of the range
        rebinned_phaii = self.phaii.rebin_time(combine_by_factor, 2, 
                                               time_range=(0.1, 0.2))
        self.assertEqual(rebinned_phaii.data.num_times, 4)
        
    def test_slice_energy(self):
        # one slice
        sliced_phaii = self.phaii.slice_energy((25.0, 750.0))
        self.assertTupleEqual(sliced_phaii.energy_range, (11.464164, 997.2431))
        self.assertEqual(sliced_phaii.num_chans, 6)
        
        # multiple slices
        sliced_phaii = self.phaii.slice_energy([(25.0, 35.0), (550.0, 750.0)])
        self.assertTupleEqual(sliced_phaii.energy_range, (11.464164, 997.2431))
        self.assertEqual(sliced_phaii.num_chans, 3)

    def test_slice_time(self):
        # one slice
        sliced_phaii = self.phaii.slice_time((0.0, 0.1))
        self.assertTupleEqual(sliced_phaii.time_range, (0.0, 0.128))
        self.assertEqual(sliced_phaii.data.num_times, 3)
        
        # multiple slices
        sliced_phaii = self.phaii.slice_time([(0.0, 0.1), (0.2, 0.3)])
        self.assertTupleEqual(sliced_phaii.time_range, (0.0, 0.32))
        self.assertEqual(sliced_phaii.data.num_times, 5)

    def test_to_lightcurve(self):
        # integrate full ranges
        lc = self.phaii.to_lightcurve()
        self.assertIsInstance(lc, TimeBins)
        self.assertTupleEqual(lc.range, (0.0, 0.320))
        
        # only create a slice of the lc
        lc = self.phaii.to_lightcurve(time_range=(0.1, 0.2))
        self.assertTupleEqual(lc.range, (0.064, 0.256))
        
        # integrate over a part of the energy range
        lc = self.phaii.to_lightcurve(energy_range=(50.0, 300.0))
        self.assertListEqual(lc.counts.tolist(), [3, 31, 29, 31, 28, 23])
        
        # integrate over the same channel range
        lc = self.phaii.to_lightcurve(channel_range=(3,5))
        self.assertListEqual(lc.counts.tolist(), [3, 31, 29, 31, 28, 23])

    def test_to_spectrum(self):
        # integrate full ranges
        spec = self.phaii.to_spectrum()
        self.assertIsInstance(spec, EnergyBins)
        self.assertTupleEqual(spec.range, (4.323754, 2000.))

        # only create a slice of the spec
        spec = self.phaii.to_spectrum(energy_range=(50.0, 300.0))
        self.assertTupleEqual(spec.range, (49.60019, 538.1436))
        
        # or channels
        spec = self.phaii.to_spectrum(channel_range=(3, 5))
        self.assertTupleEqual(spec.range, (49.60019, 538.1436))
        
        # integrate over a part of the time range
        spec = self.phaii.to_spectrum(time_range=(0.0, 0.1))
        self.assertListEqual(spec.counts.tolist(), [6, 39, 38, 27, 24, 12, 8, 8])
    
    def test_to_pha(self):
        # full time and energy range
        pha = self.phaii.to_pha()
        self.assertIsInstance(pha, Pha)
        self.assertAlmostEqual(pha.exposure, 0.3188)
        self.assertTupleEqual(pha.time_range, (0.0, 0.320))
        self.assertEqual(pha.num_chans, 8)
    
        # integrate over time slice
        pha = self.phaii.to_pha(time_ranges=(0.0, 0.1))
        self.assertAlmostEqual(pha.exposure, 0.1274)
        self.assertTupleEqual(pha.time_range, (0.0, 0.128))
        self.assertEqual(pha.num_chans, 8)
        
        # subset of the energy range
        pha = self.phaii.to_pha(energy_range=(50.0, 300.0))
        self.assertAlmostEqual(pha.exposure, 0.3188)
        self.assertTupleEqual(pha.time_range, (0.0, 0.320))
        self.assertEqual(pha.num_chans, 8)
        self.assertListEqual(pha.valid_channels.tolist(), [3, 4, 5])

        # subset of the channel range
        pha = self.phaii.to_pha(channel_range=(3, 5))
        self.assertAlmostEqual(pha.exposure, 0.3188)
        self.assertTupleEqual(pha.time_range, (0.0, 0.320))
        self.assertEqual(pha.num_chans, 8)
        self.assertListEqual(pha.valid_channels.tolist(), [3, 4, 5])

    def test_write(self):        
        with self.assertRaises(NameError):
            self.phaii.write('.')

    def test_merge(self):
        phaii1 = self.phaii.slice_time((0.0, 0.1))
        phaii2 = self.phaii.slice_time((0.15, 0.3))
        phaii = Phaii.merge([phaii1, phaii2])
        self.assertListEqual(phaii.data.tstart.tolist(), 
                             self.phaii.data.tstart.tolist())
        self.assertListEqual(phaii.data.tstop.tolist(), 
                             self.phaii.data.tstop.tolist())
        self.assertListEqual(phaii.data.emin.tolist(), 
                             self.phaii.data.emin.tolist())
        self.assertListEqual(phaii.data.emax.tolist(), 
                             self.phaii.data.emax.tolist())
        self.assertListEqual(phaii.data.exposure.tolist(), 
                             self.phaii.data.exposure.tolist())
        self.assertTupleEqual(phaii.gti.as_list()[0], 
                             self.phaii.gti.as_list()[0])
        
        # not a list of valid PHaii
        with self.assertRaises(ValueError):
            Phaii.merge([phaii1, phaii2.gti])
        
        # not same number of energy channels
        with self.assertRaises(ValueError):
            Phaii.merge([phaii1, phaii2.slice_energy((50.0, 300.0))])

    def test_no_gti(self):
        phaii = Phaii.from_data(self.phaii.data)
        self.assertTupleEqual(phaii.gti.as_list()[0], (0.0, 0.320))
    
    def test_errors(self):
        
        with self.assertRaises(TypeError):
            Phaii.from_data(self.phaii.gti)

        with self.assertRaises(TypeError):
            Phaii.from_data(self.phaii.data, gti=self.phaii.data)
    
        with self.assertRaises(ValueError):
            Phaii.from_data(self.phaii.data, trigger_time=-10.0)

        with self.assertRaises(TypeError):
            Phaii.from_data(self.phaii.data, headers=self.phaii.data)
        
        with self.assertRaises(AssertionError):
            self.phaii.slice_time((0.1, 0.0))
            

class TestPhaiiNoEbounds(unittest.TestCase):
    
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
        chan_nums = [0, 1, 2, 3, 4, 5, 6, 7]

        data = TimeChannelBins(counts, tstart, tstop, exposure, chan_nums)
        gti = Gti.from_list([(0.0000, 0.320)])
        cls.phaii = Phaii.from_data(data, gti=gti, trigger_time=356223561.133346)

    def test_data(self):
        self.assertIsInstance(self.phaii.data, TimeChannelBins)
    
    def test_ebounds(self):
        assert self.phaii.ebounds is None
    
    def test_energy_range(self):
        assert self.phaii.energy_range is None

    def test_filename(self):
        self.assertIsNone(self.phaii.filename)
    
    def test_gti(self):
        self.assertIsInstance(self.phaii.gti, Gti)
    
    def test_headers(self):
        self.assertIsNone(self.phaii.headers)

    def test_num_chans(self):
        self.assertEqual(self.phaii.num_chans, 8)
    
    def test_time_range(self):
        self.assertTupleEqual(self.phaii.time_range, (0.0, 0.320))
    
    def test_trigtime(self):
        self.assertEqual(self.phaii.trigtime, 356223561.133346)

    def test_get_exposure(self):
        # total exposure
        self.assertAlmostEqual(self.phaii.get_exposure(), 0.3188)
        
        # exposure of a time slice
        self.assertAlmostEqual(self.phaii.get_exposure(time_ranges=(0.0, 0.1)),
                               0.1274)

        # exposure of multiple time slices
        self.assertAlmostEqual(self.phaii.get_exposure(time_ranges=[(0.0, 0.1),
                                                                    (0.2, 0.3)]),
                               0.255)
    
    def test_rebin_energy(self):
        # full range
        rebinned_phaii = self.phaii.rebin_energy(combine_by_factor, 2)
        self.assertEqual(rebinned_phaii.num_chans, 4)
        
        # part of the range
        rebinned_phaii = self.phaii.rebin_energy(combine_by_factor, 2, 
                                                 energy_range=(1, 6))
        self.assertEqual(rebinned_phaii.num_chans, 5)

    def test_rebin_time(self):
        # full range
        rebinned_phaii = self.phaii.rebin_time(combine_by_factor, 2)
        self.assertEqual(rebinned_phaii.data.num_times, 3)

        # part of the range
        rebinned_phaii = self.phaii.rebin_time(combine_by_factor, 2, 
                                               time_range=(0.1, 0.2))
        self.assertEqual(rebinned_phaii.data.num_times, 4)
        
    def test_slice_energy(self):
        # one slice
        sliced_phaii = self.phaii.slice_energy((1, 6))
        assert sliced_phaii.energy_range is None
        self.assertEqual(sliced_phaii.num_chans, 6)
        
        # multiple slices
        sliced_phaii = self.phaii.slice_energy([(1, 2), (6, 6)])
        assert sliced_phaii.energy_range is None
        self.assertEqual(sliced_phaii.num_chans, 3)

    def test_slice_time(self):
        # one slice
        sliced_phaii = self.phaii.slice_time((0.0, 0.1))
        self.assertTupleEqual(sliced_phaii.time_range, (0.0, 0.128))
        self.assertEqual(sliced_phaii.data.num_times, 3)
        
        # multiple slices
        sliced_phaii = self.phaii.slice_time([(0.0, 0.1), (0.2, 0.3)])
        self.assertTupleEqual(sliced_phaii.time_range, (0.0, 0.32))
        self.assertEqual(sliced_phaii.data.num_times, 5)

    def test_to_lightcurve(self):
        # integrate full ranges
        lc = self.phaii.to_lightcurve()
        self.assertIsInstance(lc, TimeBins)
        self.assertTupleEqual(lc.range, (0.0, 0.320))
        
        # only create a slice of the lc
        lc = self.phaii.to_lightcurve(time_range=(0.1, 0.2))
        self.assertTupleEqual(lc.range, (0.064, 0.256))
        
        # this should ignore the energy range and integrate over full range
        # because we have data with no energy calibration
        lc = self.phaii.to_lightcurve(energy_range=(50.0, 300.0))
        self.assertListEqual(lc.counts.tolist(), [5, 66, 91, 82, 75, 76])
        
        # integrate over the same channel range
        lc = self.phaii.to_lightcurve(channel_range=(3,5))
        self.assertListEqual(lc.counts.tolist(), [3, 31, 29, 31, 28, 23])

    def test_to_spectrum(self):
        # integrate full ranges
        spec = self.phaii.to_spectrum()
        self.assertIsInstance(spec, ChannelBins)
        self.assertTupleEqual(spec.range, (0, 7))

        # this should ignore the energy range and use the full range
        # because we have data with no energy calibration
        spec = self.phaii.to_spectrum(energy_range=(50.0, 300.0))
        self.assertTupleEqual(spec.range, (0, 7))
        
        # or channels
        spec = self.phaii.to_spectrum(channel_range=(3, 5))
        self.assertTupleEqual(spec.range, (3, 5))
        
        # integrate over a part of the time range
        spec = self.phaii.to_spectrum(time_range=(0.0, 0.1))
        self.assertListEqual(spec.counts.tolist(), [6, 39, 38, 27, 24, 12, 8, 8])
    
    def test_to_pha(self):
        
        # should raise an exception because we can't create a PHA with no 
        # energy calibration
        with self.assertRaises(RuntimeError):
            self.phaii.to_pha()

    def test_write(self):        
        with self.assertRaises(NameError):
            self.phaii.write('.')

    def test_merge(self):
        phaii1 = self.phaii.slice_time((0.0, 0.1))
        phaii2 = self.phaii.slice_time((0.15, 0.3))
        phaii = Phaii.merge([phaii1, phaii2])
        self.assertListEqual(phaii.data.tstart.tolist(), 
                             self.phaii.data.tstart.tolist())
        self.assertListEqual(phaii.data.tstop.tolist(), 
                             self.phaii.data.tstop.tolist())
        self.assertListEqual(phaii.data.chan_nums.tolist(), 
                             self.phaii.data.chan_nums.tolist())
        self.assertListEqual(phaii.data.exposure.tolist(), 
                             self.phaii.data.exposure.tolist())
        self.assertTupleEqual(phaii.gti.as_list()[0], 
                             self.phaii.gti.as_list()[0])
        
        # not a list of valid PHaii
        with self.assertRaises(ValueError):
            Phaii.merge([phaii1, phaii2.gti])
        
        # not same number of energy channels
        with self.assertRaises(ValueError):
            Phaii.merge([phaii1, phaii2.slice_energy((3, 4))])

    def test_no_gti(self):
        phaii = Phaii.from_data(self.phaii.data)
        self.assertTupleEqual(phaii.gti.as_list()[0], (0.0, 0.320))
    
    def test_errors(self):
        
        with self.assertRaises(TypeError):
            Phaii.from_data(self.phaii.gti)

        with self.assertRaises(TypeError):
            Phaii.from_data(self.phaii.data, gti=self.phaii.data)
    
        with self.assertRaises(ValueError):
            Phaii.from_data(self.phaii.data, trigger_time=-10.0)

        with self.assertRaises(TypeError):
            Phaii.from_data(self.phaii.data, headers=self.phaii.data)
        
        with self.assertRaises(AssertionError):
            self.phaii.slice_time((0.1, 0.0))
            
                
if __name__ == '__main__':
    unittest.main()


