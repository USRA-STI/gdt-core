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
import numpy as np
from gdt.core.data_primitives import EventList, Gti, Ebounds, EnergyBins, ChannelBins
from gdt.core.headers import FileHeaders
from gdt.core.phaii import Phaii
from gdt.core.pha import Pha
from gdt.core.tte import PhotonList
from gdt.core.binning.binned import combine_by_factor
from gdt.core.binning.unbinned import bin_by_time

class TestPhotonList(unittest.TestCase):
    
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
        
    def test_data(self):
        self.assertIsInstance(self.tte.data, EventList)
    
    def test_ebounds(self):
        self.assertIsInstance(self.tte.ebounds, Ebounds)
    
    def test_energy_range(self):
        self.assertTupleEqual(self.tte.energy_range, (10.0, 640.))

    def test_event_deadtime(self):
        self.assertEqual(self.tte.event_deadtime, 0.001)

    def test_filename(self):
        self.assertIsNone(self.tte.filename)
    
    def test_gti(self):
        self.assertIsInstance(self.tte.gti, Gti)
    
    def test_headers(self):
        self.assertIsNone(self.tte.headers)

    def test_num_chans(self):
        self.assertEqual(self.tte.num_chans, 6)

    def test_overflow_deadtime(self):
        self.assertEqual(self.tte.overflow_deadtime, 0.1)
    
    def test_time_range(self):
        self.assertTupleEqual(self.tte.time_range, (0.706, 10.61))
    
    def test_trigtime(self):
        self.assertEqual(self.tte.trigtime, 356223561.133346)

    def test_get_exposure(self):
        # total exposure
        self.assertAlmostEqual(self.tte.get_exposure(), 9.795)
        
        # exposure of a time slice
        self.assertAlmostEqual(self.tte.get_exposure(time_ranges=(0.0, 5.0)),
                               4.996)

        # exposure of multiple time slices
        self.assertAlmostEqual(self.tte.get_exposure(time_ranges=[(0.0, 2.0),
                                                                    (5.0, 7.0)]),
                               3.898)
    
    def test_rebin_energy(self):
        # full range
        rebinned_tte = self.tte.rebin_energy(combine_by_factor, 2)
        self.assertEqual(rebinned_tte.num_chans, 3)        
        
    def test_slice_energy(self):
        # one slice
        sliced_tte = self.tte.slice_energy((50.0, 250.0))
        self.assertTupleEqual(sliced_tte.energy_range, (40.0, 320.0))
        self.assertEqual(sliced_tte.num_chans, 6)
        
        # multiple slices
        sliced_tte = self.tte.slice_energy([(25.0, 75.0), (250.0, 500.0)])
        self.assertTupleEqual(sliced_tte.energy_range, (20.0, 640.0))
        self.assertEqual(sliced_tte.num_chans, 6)

    def test_slice_time(self):
        # one slice
        sliced_tte = self.tte.slice_time((0.0, 5.0))
        self.assertTupleEqual(sliced_tte.time_range, (0.706, 3.512))
        self.assertEqual(sliced_tte.data.size, 4)
        
        # multiple slices
        sliced_tte = self.tte.slice_time([(0.0, 2.0), (5.0, 7.0)])
        self.assertTupleEqual(sliced_tte.time_range, (0.706, 5.540))
        self.assertEqual(sliced_tte.data.size, 3)

    def test_to_spectrum(self):
        # integrate full ranges
        spec = self.tte.to_spectrum()
        self.assertIsInstance(spec, EnergyBins)
        self.assertTupleEqual(spec.range, (10.0, 640.))

        # only create a slice of the spec
        spec = self.tte.to_spectrum(energy_range=(50.0, 300.0))
        self.assertListEqual(spec.counts.tolist(), [0, 0, 1, 0, 3, 0])
        
        # or channels
        spec = self.tte.to_spectrum(channel_range=(3, 5))
        self.assertListEqual(spec.counts.tolist(), [0, 0, 0, 0, 3, 1])
        
        # integrate over a part of the time range
        spec = self.tte.to_spectrum(time_range=(0.0, 5.0))
        self.assertListEqual(spec.counts.tolist(), [1, 1, 0, 0, 2, 0])
    
    def test_to_pha(self):
        
        # full time and energy range
        pha = self.tte.to_pha()
        self.assertIsInstance(pha, Pha)
        self.assertAlmostEqual(pha.exposure, 9.795)
        self.assertTupleEqual(pha.time_range, (0.706, 10.61))
        self.assertEqual(pha.num_chans, 6)
        
        # integrate over time slice
        pha = self.tte.to_pha(time_ranges=(0.0, 5.0))
        self.assertAlmostEqual(pha.exposure, 2.802)
        self.assertTupleEqual(pha.time_range, (0.706, 3.512))
        self.assertEqual(pha.num_chans, 6)
        
        # subset of the energy range
        pha = self.tte.to_pha(energy_range=(50.0, 300.0))
        self.assertAlmostEqual(pha.exposure, 9.74)
        self.assertTupleEqual(pha.time_range, (0.706, 10.61))
        self.assertEqual(pha.num_chans, 6)
        self.assertListEqual(pha.valid_channels.tolist(), [2, 3, 4])
        
        # subset of the channel range
        pha = self.tte.to_pha(channel_range=(3, 5))
        self.assertAlmostEqual(pha.exposure, 8.793)
        self.assertTupleEqual(pha.time_range, (0.706, 10.61))
        self.assertEqual(pha.num_chans, 6)
        self.assertListEqual(pha.valid_channels.tolist(), [3, 4, 5])
    
    def test_to_phaii(self):
        # full range
        phaii = self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii)
        self.assertTupleEqual(phaii.energy_range, self.tte.energy_range)
        self.assertTupleEqual(phaii.gti.as_list()[0], self.tte.gti.as_list()[0])
        self.assertEqual(phaii.num_chans, 6)
        self.assertTupleEqual(phaii.time_range, (0.706, 10.706))
        self.assertEqual(phaii.trigtime, self.tte.trigtime)
        self.assertEqual(phaii.data.num_times, 10)
        
        # time range
        phaii = self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, 
                                  time_range=(0.0, 5.0))
        self.assertTupleEqual(phaii.energy_range, self.tte.energy_range)
        self.assertTupleEqual(phaii.gti.as_list()[0], (0.706, 3.512))
        self.assertEqual(phaii.num_chans, 6)
        self.assertTupleEqual(phaii.time_range, (0.706, 3.706))
        self.assertEqual(phaii.trigtime, self.tte.trigtime)
        self.assertEqual(phaii.data.num_times, 3)
        
        # energy range
        phaii = self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, 
                                  energy_range=(10.0, 50.0))
        self.assertTupleEqual(phaii.energy_range, (10.0, 80.0))
        self.assertTupleEqual(phaii.gti.as_list()[0], self.tte.gti.as_list()[0])
        self.assertEqual(phaii.num_chans, 3)
        self.assertEqual(phaii.trigtime, self.tte.trigtime)
        self.assertEqual(phaii.data.num_times, 9)

        # channel range
        phaii = self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, 
                                  channel_range=(0, 2))
        self.assertTupleEqual(phaii.energy_range, (10.0, 80.0))
        self.assertTupleEqual(phaii.gti.as_list()[0], self.tte.gti.as_list()[0])
        self.assertEqual(phaii.num_chans, 3)
        self.assertEqual(phaii.trigtime, self.tte.trigtime)
        self.assertEqual(phaii.data.num_times, 9)
        
        # bad time range
        with self.assertRaises(AssertionError):
            self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, 
                              time_range=(1.0, 0.0))

        # bad energy range
        with self.assertRaises(AssertionError):
            self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, 
                              energy_range=(1.0, 0.0))
       
        # bad channel range
        with self.assertRaises(AssertionError):
            self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, 
                              channel_range=(1.0, 0.0))
        
        # bad Phaii class
        with self.assertRaises(TypeError):
            self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Pha)
    
    def test_write(self):
        with self.assertRaises(NameError):
            self.tte.write('.')

    def test_merge(self):
        tte1 = self.tte.slice_time((0.0, 5.0))
        tte2 = self.tte.slice_time((5.5, 11.0))
        tte = PhotonList.merge([tte1, tte2])
        self.assertListEqual(tte.data.times.tolist(), 
                             self.tte.data.times.tolist())
        self.assertListEqual(tte.data.channels.tolist(), 
                             self.tte.data.channels.tolist())
        self.assertTupleEqual(tte.gti.as_list()[0], (0.706, 3.512))
        self.assertTupleEqual(tte.gti.as_list()[1], (5.54, 10.61))
        
        # not a list of valid PhotonLists
        with self.assertRaises(ValueError):
            PhotonList.merge([tte1, tte2.gti])
        
        # not a valid header index
        with self.assertRaises(ValueError):
            PhotonList.merge([tte1, tte2], primary=2)
        
        eb = Ebounds.from_list(tte2.ebounds.as_list()[:-1])
        el = EventList(tte2.data.times, tte2.data.channels, ebounds=eb)
        tte3 = PhotonList.from_data(el, gti=tte2.gti, 
                                    trigger_time=tte2.trigtime)
        with self.assertRaises(ValueError):
            PhotonList.merge([tte1, tte3])
        
    def test_no_gti(self):
        tte = PhotonList.from_data(self.tte.data)
        self.assertTupleEqual(tte.gti.as_list()[0], (0.706, 10.61))
    
    def test_errors(self):
        
        with self.assertRaises(TypeError):
            PhotonList.from_data(self.tte.gti)

        with self.assertRaises(TypeError):
            PhotonList.from_data(self.tte.data, gti=self.tte.data)
    
        with self.assertRaises(ValueError):
            PhotonList.from_data(self.tte.data, trigger_time=-10.0)

        with self.assertRaises(TypeError):
            PhotonList.from_data(self.tte.data, headers=self.tte.data)

        with self.assertRaises(ValueError):
            PhotonList.from_data(self.tte.data, event_deadtime=-10.0)

        with self.assertRaises(ValueError):
            PhotonList.from_data(self.tte.data, overflow_deadtime=-10.0)

        with self.assertRaises(TypeError):
            PhotonList.from_data(self.tte.data, event_deadtime='')

        with self.assertRaises(TypeError):
            PhotonList.from_data(self.tte.data, overflow_deadtime='')


class TestPhotonListNoEbounds(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        times = [0.706, 1.640, 3.185, 3.512, 5.540, 
                 7.790, 9.602, 9.726, 10.45, 10.61]
        chans = [4, 1, 0, 4, 5, 0, 4, 0, 2, 0]
        data = EventList(times=times, channels=chans)
        gti = Gti.from_list([(0.0000, 10.70)])
        cls.tte = PhotonList.from_data(data, gti=gti, 
                                       trigger_time=356223561.133346,
                                       event_deadtime=0.001, 
                                       overflow_deadtime=0.1)
        
    def test_data(self):
        self.assertIsInstance(self.tte.data, EventList)
    
    def test_ebounds(self):
        assert self.tte.ebounds is None
    
    def test_energy_range(self):
        assert self.tte.energy_range is None

    def test_event_deadtime(self):
        self.assertEqual(self.tte.event_deadtime, 0.001)

    def test_filename(self):
        self.assertIsNone(self.tte.filename)
    
    def test_gti(self):
        self.assertIsInstance(self.tte.gti, Gti)
    
    def test_headers(self):
        self.assertIsNone(self.tte.headers)

    def test_num_chans(self):
        self.assertEqual(self.tte.num_chans, 6)

    def test_overflow_deadtime(self):
        self.assertEqual(self.tte.overflow_deadtime, 0.1)
    
    def test_time_range(self):
        self.assertTupleEqual(self.tte.time_range, (0.706, 10.61))
    
    def test_trigtime(self):
        self.assertEqual(self.tte.trigtime, 356223561.133346)

    def test_get_exposure(self):
        # total exposure
        self.assertAlmostEqual(self.tte.get_exposure(), 9.795)
        
        # exposure of a time slice
        self.assertAlmostEqual(self.tte.get_exposure(time_ranges=(0.0, 5.0)),
                               4.996)

        # exposure of multiple time slices
        self.assertAlmostEqual(self.tte.get_exposure(time_ranges=[(0.0, 2.0),
                                                                    (5.0, 7.0)]),
                               3.898)
    
    def test_rebin_energy(self):
        # full range
        rebinned_tte = self.tte.rebin_energy(combine_by_factor, 2)
        self.assertEqual(rebinned_tte.num_chans, 3)        
        
    def test_slice_energy(self):
        # one slice
        sliced_tte = self.tte.slice_energy((2, 4))
        assert sliced_tte.energy_range is None
        self.assertEqual(sliced_tte.num_chans, 3)
        
        # multiple slices
        sliced_tte = self.tte.slice_energy([(1, 2), (4, 5)])
        assert sliced_tte.energy_range is None
        self.assertEqual(sliced_tte.num_chans, 5)

    def test_slice_time(self):
        # one slice
        sliced_tte = self.tte.slice_time((0.0, 5.0))
        self.assertTupleEqual(sliced_tte.time_range, (0.706, 3.512))
        self.assertEqual(sliced_tte.data.size, 4)
        
        # multiple slices
        sliced_tte = self.tte.slice_time([(0.0, 2.0), (5.0, 7.0)])
        self.assertTupleEqual(sliced_tte.time_range, (0.706, 5.540))
        self.assertEqual(sliced_tte.data.size, 3)

    def test_to_spectrum(self):
        # integrate full ranges
        spec = self.tte.to_spectrum()
        self.assertIsInstance(spec, ChannelBins)
        self.assertTupleEqual(spec.range, (0, 5))
        
        # this should be the whole range because we have no ebounds
        spec = self.tte.to_spectrum(energy_range=(50.0, 300.0))
        self.assertListEqual(spec.counts.tolist(), [4, 1, 1, 0, 3, 1])
        
        # or channels
        spec = self.tte.to_spectrum(channel_range=(3, 5))
        self.assertListEqual(spec.counts.tolist(), [0, 0, 0, 0, 3, 1])
        
        # integrate over a part of the time range
        spec = self.tte.to_spectrum(time_range=(0.0, 5.0))
        self.assertListEqual(spec.counts.tolist(), [1, 1, 0, 0, 2])
    
    def test_to_pha(self):
        
        # this should raise an exception
        with self.assertRaises(RuntimeError):
            self.tte.to_pha()
        
    def test_to_phaii(self):
        # full range
        phaii = self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii)
        assert phaii.energy_range is None
        self.assertTupleEqual(phaii.gti.as_list()[0], self.tte.gti.as_list()[0])
        self.assertEqual(phaii.num_chans, 6)
        self.assertTupleEqual(phaii.time_range, (0.706, 10.706))
        self.assertEqual(phaii.trigtime, self.tte.trigtime)
        self.assertEqual(phaii.data.num_times, 10)
        
        # time range
        phaii = self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, 
                                  time_range=(0.0, 5.0))
        assert phaii.energy_range is None
        self.assertTupleEqual(phaii.gti.as_list()[0], (0.706, 3.512))
        self.assertEqual(phaii.num_chans, 5)
        self.assertTupleEqual(phaii.time_range, (0.706, 3.706))
        self.assertEqual(phaii.trigtime, self.tte.trigtime)
        self.assertEqual(phaii.data.num_times, 3)
        
        # energy range (this should be same as full range because no ebounds)
        phaii = self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, 
                                  energy_range=(10.0, 50.0))
        assert phaii.energy_range is None
        self.assertTupleEqual(phaii.gti.as_list()[0], self.tte.gti.as_list()[0])
        self.assertEqual(phaii.num_chans, 6)
        self.assertEqual(phaii.trigtime, self.tte.trigtime)
        self.assertEqual(phaii.data.num_times, 10)

        # channel range
        phaii = self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, 
                                  channel_range=(0, 2))
        assert phaii.energy_range is None
        self.assertTupleEqual(phaii.gti.as_list()[0], self.tte.gti.as_list()[0])
        self.assertEqual(phaii.num_chans, 3)
        self.assertEqual(phaii.trigtime, self.tte.trigtime)
        self.assertEqual(phaii.data.num_times, 9)
        
        # bad time range
        with self.assertRaises(AssertionError):
            self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, 
                              time_range=(1.0, 0.0))
       
        # bad channel range
        with self.assertRaises(AssertionError):
            self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, 
                              channel_range=(1.0, 0.0))
        
        # bad Phaii class
        with self.assertRaises(TypeError):
            self.tte.to_phaii(bin_by_time, 1.0, phaii_class=Pha)
    
    def test_write(self):
        with self.assertRaises(NameError):
            self.tte.write('.')

    def test_merge(self):
        tte1 = self.tte.slice_time((0.0, 5.0))
        tte2 = self.tte.slice_time((6.0, 10.0))
        tte = PhotonList.merge([tte1, tte2])
        self.assertListEqual(tte.data.times.tolist(), 
                             [0.706, 1.640, 3.185, 3.512, 7.790, 9.602, 9.726])
        self.assertListEqual(tte.data.channels.tolist(), 
                             [4, 1, 0, 4, 0, 4, 0])
        self.assertTupleEqual(tte.gti.as_list()[0], (0.706, 3.512))
        self.assertTupleEqual(tte.gti.as_list()[1], (7.79, 9.726))
        
        # not a list of valid PhotonLists
        with self.assertRaises(ValueError):
            PhotonList.merge([tte1, tte2.gti])
        
        # not a valid header index
        with self.assertRaises(ValueError):
            PhotonList.merge([tte1, tte2], primary=2)
                
    def test_no_gti(self):
        tte = PhotonList.from_data(self.tte.data)
        self.assertTupleEqual(tte.gti.as_list()[0], (0.706, 10.61))
    
    def test_errors(self):
        
        with self.assertRaises(TypeError):
            PhotonList.from_data(self.tte.gti)

        with self.assertRaises(TypeError):
            PhotonList.from_data(self.tte.data, gti=self.tte.data)
    
        with self.assertRaises(ValueError):
            PhotonList.from_data(self.tte.data, trigger_time=-10.0)

        with self.assertRaises(TypeError):
            PhotonList.from_data(self.tte.data, headers=self.tte.data)

        with self.assertRaises(ValueError):
            PhotonList.from_data(self.tte.data, event_deadtime=-10.0)

        with self.assertRaises(ValueError):
            PhotonList.from_data(self.tte.data, overflow_deadtime=-10.0)

        with self.assertRaises(TypeError):
            PhotonList.from_data(self.tte.data, event_deadtime='')

        with self.assertRaises(TypeError):
            PhotonList.from_data(self.tte.data, overflow_deadtime='')
                    
                
if __name__ == '__main__':
    unittest.main()


