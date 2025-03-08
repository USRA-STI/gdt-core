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
from gdt.core.data_primitives import EnergyBins, Gti, Ebounds
from gdt.core.headers import FileHeaders
from gdt.core.pha import Pha, Bak
from gdt.core.binning.binned import combine_by_factor
from gdt.core.background.primitives import BackgroundSpectrum

this_dir = os.path.dirname(__file__)


class TestPha(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        counts = [119, 71, 34, 30, 21, 6, 2, 19]
        emin = [4.323754, 11.464164, 26.22962, 49.60019, 101.016815,
                290.46063, 538.1436, 997.2431]
        emax = [11.464164, 26.22962, 49.60019, 101.016815, 290.46063,
                538.1436, 997.2431, 2000.]
        exposure = 0.25459924
         
        data = EnergyBins(counts, emin, emax, exposure)
        gti = Gti.from_list([(-899.0864419937134, -898.8306360244751)])
        cls.pha = Pha.from_data(data, gti=gti, trigger_time=356223561.133346)

    @classmethod
    def tearDownClass(cls):
        try:
            os.remove(os.path.join(this_dir, 'test_pha.pha'))
        except:
            pass
        cls.pha.close()

    def test_channel_mask(self):
        self.assertListEqual(self.pha.channel_mask.tolist(), [True]*8)
        
        chan_mask = [False, False, False, False, True, True, True, True]
        new_pha = Pha.from_data(self.pha.data, gti=self.pha.gti, 
                                trigger_time=self.pha.trigtime,
                                channel_mask=chan_mask)
        self.assertListEqual(new_pha.channel_mask.tolist(), chan_mask)

    def test_data(self):
        self.assertIsInstance(self.pha.data, EnergyBins)
    
    def test_ebounds(self):
        self.assertIsInstance(self.pha.ebounds, Ebounds)
    
    def test_energy_range(self):
        self.assertTupleEqual(self.pha.energy_range, (4.323754, 2000.))

    def test_exposure(self):
        self.assertEqual(self.pha.exposure, 0.25459924)
    
    def test_filename(self):
        self.assertIsNone(self.pha.filename)
    
    def test_gti(self):
        self.assertIsInstance(self.pha.gti, Gti)
    
    def test_headers(self):
        self.assertIsInstance(self.pha.headers, FileHeaders)
        self.assertEqual(self.pha.headers['PRIMARY']['TSTART'], -899.0864419937134)
        self.assertEqual(self.pha.headers['PRIMARY']['TSTOP'], -898.8306360244751)
        self.assertEqual(self.pha.headers['PRIMARY']['TRIGTIME'], 356223561.133346)
        self.assertEqual(self.pha.headers['EBOUNDS']['DETCHANS'], 8)
        self.assertEqual(self.pha.headers['SPECTRUM']['DETCHANS'], 8)
        self.assertEqual(self.pha.headers['SPECTRUM']['EXPOSURE'], 0.25459924)

    def test_num_chans(self):
        self.assertEqual(self.pha.num_chans, 8)
    
    def test_tcent(self):
        self.assertAlmostEqual(self.pha.tcent, -898.958539)
    
    def test_time_range(self):
        self.assertTupleEqual(self.pha.time_range, (-899.0864419937134, 
                                                    -898.8306360244751))
    
    def test_trigtime(self):
        self.assertEqual(self.pha.trigtime, 356223561.133346)

    def test_valid_channels(self):
        self.assertListEqual(self.pha.valid_channels.tolist(), [0,1,2,3,4,5,6,7])

        chan_mask = [False, False, False, False, True, True, True, True]
        new_pha = Pha.from_data(self.pha.data, gti=self.pha.gti, 
                                trigger_time=self.pha.trigtime,
                                channel_mask=chan_mask)
        self.assertListEqual(new_pha.valid_channels.tolist(), [4, 5, 6, 7])

    def test_rebin_energy(self):
        # full range
        rebinned_pha = self.pha.rebin_energy(combine_by_factor, 2)
        self.assertEqual(rebinned_pha.num_chans, 4)
        self.assertEqual(rebinned_pha.headers['EBOUNDS']['DETCHANS'], 4)
        self.assertEqual(rebinned_pha.headers['SPECTRUM']['DETCHANS'], 4)
        
        # part of the range
        rebinned_pha = self.pha.rebin_energy(combine_by_factor, 2, 
                                             energy_range=(25.0, 750.0))
        self.assertEqual(rebinned_pha.num_chans, 5)
        self.assertEqual(rebinned_pha.headers['EBOUNDS']['DETCHANS'], 5)
        self.assertEqual(rebinned_pha.headers['SPECTRUM']['DETCHANS'], 5)

    def test_slice_energy(self):
        # one slice
        sliced_pha = self.pha.slice_energy((25.0, 750.0))
        self.assertTupleEqual(sliced_pha.energy_range, (11.464164, 997.2431))
        self.assertEqual(sliced_pha.num_chans, 6)
        self.assertEqual(sliced_pha.headers['EBOUNDS']['DETCHANS'], 6)
        self.assertEqual(sliced_pha.headers['SPECTRUM']['DETCHANS'], 6)
        
        # multiple slices
        sliced_pha = self.pha.slice_energy([(25.0, 35.0), (550.0, 750.0)])
        self.assertTupleEqual(sliced_pha.energy_range, (11.464164, 997.2431))
        self.assertEqual(sliced_pha.num_chans, 3)
        self.assertEqual(sliced_pha.headers['EBOUNDS']['DETCHANS'], 3)
        self.assertEqual(sliced_pha.headers['SPECTRUM']['DETCHANS'], 3)

    def test_write_read(self):
        self.pha.write(this_dir, 'test_pha.pha')
        pha2 = Pha.open(os.path.join(this_dir, 'test_pha.pha'))

        # test data
        self.assertListEqual(pha2.data.counts.tolist(), 
                             self.pha.data.counts.tolist())
        for i in range(pha2.num_chans):
            self.assertAlmostEqual(pha2.data.lo_edges[i],
                                   self.pha.data.lo_edges[i], places=4)
        for i in range(pha2.num_chans):
            self.assertAlmostEqual(pha2.data.hi_edges[i],
                                   self.pha.data.hi_edges[i], places=4)
        # test ebounds
        for i in range(pha2.num_chans):
            self.assertAlmostEqual(pha2.ebounds.low_edges()[i],
                                   self.pha.ebounds.low_edges()[i], places=4)
        for i in range(pha2.num_chans):
            self.assertAlmostEqual(pha2.ebounds.high_edges()[i],
                                   self.pha.ebounds.high_edges()[i], places=4)

        # test gti
        self.assertListEqual(pha2.gti.low_edges(), self.pha.gti.low_edges())
        self.assertListEqual(pha2.gti.high_edges(), self.pha.gti.high_edges())

        # test attributes
        self.assertEqual(pha2.filename, 'test_pha.pha')
        self.assertEqual(pha2.num_chans, self.pha.num_chans)
        self.assertEqual(pha2.trigtime, self.pha.trigtime)
        
        # test headers
        self.assertEqual(pha2.headers['PRIMARY']['TSTART'], -899.0864419937134)
        self.assertEqual(pha2.headers['PRIMARY']['TSTOP'], -898.8306360244751)
        self.assertEqual(pha2.headers['PRIMARY']['TRIGTIME'], 356223561.133346)
        self.assertEqual(pha2.headers['EBOUNDS']['DETCHANS'], 8)
        self.assertEqual(pha2.headers['SPECTRUM']['DETCHANS'], 8)
        self.assertEqual(pha2.headers['SPECTRUM']['EXPOSURE'], 0.25459924)
        
        pha2.write(this_dir, overwrite=True)
        pha2.close()

    def test_no_gti(self):
        pha = Pha.from_data(self.pha.data)
        self.assertTupleEqual(pha.time_range, (0.0, 0.25459924))
    
    def test_errors(self):
        
        with self.assertRaises(TypeError):
            Pha.from_data(self.pha.gti)

        with self.assertRaises(TypeError):
            Pha.from_data(self.pha.data, gti=self.pha.data)
    
        with self.assertRaises(ValueError):
            Pha.from_data(self.pha.data, trigger_time=-10.0)

        with self.assertRaises(TypeError):
            Pha.from_data(self.pha.data, headers=self.pha.data)

        with self.assertRaises(TypeError):
            Pha.from_data(self.pha.data, channel_mask=self.pha.data)

        with self.assertRaises(ValueError):
            Pha.from_data(self.pha.data, channel_mask=self.pha.channel_mask[1:])
        
        with self.assertRaises(NameError):
            pha = Pha.from_data(self.pha.data)
            pha.write(this_dir)


class TestPhaChannelMask(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        counts = [119, 71, 34, 30, 21, 6, 2, 19]
        emin = [4.323754, 11.464164, 26.22962, 49.60019, 101.016815,
                290.46063, 538.1436, 997.2431]
        emax = [11.464164, 26.22962, 49.60019, 101.016815, 290.46063,
                538.1436, 997.2431, 2000.]
        exposure = 0.25459924
         
        data = EnergyBins(counts, emin, emax, exposure)
        gti = Gti.from_list([(-899.0864419937134, -898.8306360244751)])
        channel_mask = [False, True, True, True, True, True, False, False]
        cls.pha = Pha.from_data(data, gti=gti, trigger_time=356223561.133346,
                                channel_mask=channel_mask)

    def test_channel_mask(self):
        assert self.pha.channel_mask.tolist() == [False, True, True, True, 
                                                  True, True, False, False]
    
    def test_energy_range(self):
        assert self.pha.energy_range == (11.464164, 538.1436)

    def test_valid_channels(self):
        assert self.pha.valid_channels.tolist() == [1, 2, 3, 4, 5]

    def test_errors(self):
        with self.assertRaises(TypeError):
            Pha.from_data(self.pha.data, self.pha.gti, channel_mask=True)
        
        with self.assertRaises(ValueError):
            Pha.from_data(self.pha.data, self.pha.gti, channel_mask=('test',))

        with self.assertRaises(ValueError):
            Pha.from_data(self.pha.data, self.pha.gti, channel_mask=[True, False])
            


class TestPhaValidChannels(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        counts = [119, 71, 34, 30, 21, 6, 2, 19]
        emin = [4.323754, 11.464164, 26.22962, 49.60019, 101.016815,
                290.46063, 538.1436, 997.2431]
        emax = [11.464164, 26.22962, 49.60019, 101.016815, 290.46063,
                538.1436, 997.2431, 2000.]
        exposure = 0.25459924
         
        data = EnergyBins(counts, emin, emax, exposure)
        gti = Gti.from_list([(-899.0864419937134, -898.8306360244751)])
        valid_channels = [1, 2, 4, 5]
        cls.pha = Pha.from_data(data, gti=gti, trigger_time=356223561.133346,
                                valid_channels=valid_channels)

    def test_channel_mask(self):
        assert self.pha.channel_mask.tolist() == [False, True, True, False, 
                                                  True, True, False, False]
    
    def test_energy_range(self):
        assert self.pha.energy_range == (11.464164, 538.1436)

    def test_valid_channels(self):
        assert self.pha.valid_channels.tolist() == [1, 2, 4, 5]

    def test_errors(self):
        with self.assertRaises(TypeError):
            Pha.from_data(self.pha.data, self.pha.gti, valid_channels=0)
        
        with self.assertRaises(ValueError):
            Pha.from_data(self.pha.data, self.pha.gti, valid_channels=('test',))
            

class TestPhaZeroCounts(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        counts = [0, 71, 0, 30, 0, 6, 0, 0]
        emin = [4.323754, 11.464164, 26.22962, 49.60019, 101.016815,
                290.46063, 538.1436, 997.2431]
        emax = [11.464164, 26.22962, 49.60019, 101.016815, 290.46063,
                538.1436, 997.2431, 2000.]
        exposure = 0.25459924
         
        data = EnergyBins(counts, emin, emax, exposure)
        gti = Gti.from_list([(-899.0864419937134, -898.8306360244751)])
        cls.pha = Pha.from_data(data, gti=gti, trigger_time=356223561.133346)

    def test_channel_mask(self):
        assert self.pha.channel_mask.tolist() == [False, True, False, True, 
                                                  False, True, False, False]
    
    def test_energy_range(self):
        assert self.pha.energy_range == (11.464164, 538.1436)

    def test_valid_channels(self):
        assert self.pha.valid_channels.tolist() == [1, 3, 5]


class TestBak(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        rates = [467.40, 278.87, 133.54, 117.83, 82.48, 23.57, 7.86, 75.63]
        uncert = [23.37, 13.94, 6.68, 5.89, 4.12, 1.18, 0.39, 3.73]
        emin = [4.323754, 11.464164, 26.22962, 49.60019, 101.016815,
                290.46063, 538.1436, 997.2431]
        emax = [11.464164, 26.22962, 49.60019, 101.016815, 290.46063,
                538.1436, 997.2431, 2000.]
        exposure = 2.021
        
        data = BackgroundSpectrum(rates, uncert, emin, emax, exposure)
        gti = Gti.from_list([(-899.0864419937134, -897.0384419937134)])
        cls.bak = Bak.from_data(data, gti=gti, trigger_time=356223561.133346)

    @classmethod
    def tearDownClass(cls):
        try:
            os.remove(os.path.join(this_dir, 'test_bak.bak'))
        except:
            pass
        cls.bak.close()

    def test_exposure(self):
        self.assertEqual(self.bak.exposure, 2.021)

    def test_headers(self):
        self.assertIsInstance(self.bak.headers, FileHeaders)
        self.assertEqual(self.bak.headers['PRIMARY']['TSTART'], -899.0864419937134)
        self.assertEqual(self.bak.headers['PRIMARY']['TSTOP'], -897.0384419937134)
        self.assertEqual(self.bak.headers['PRIMARY']['TRIGTIME'], 356223561.133346)
        self.assertEqual(self.bak.headers['EBOUNDS']['DETCHANS'], 8)
        self.assertEqual(self.bak.headers['SPECTRUM']['DETCHANS'], 8)
        self.assertEqual(self.bak.headers['SPECTRUM']['EXPOSURE'], 2.021)

    def test_write_read(self):
        self.bak.write(this_dir, 'test_bak.bak')
        bak2 = Bak.open(os.path.join(this_dir, 'test_bak.bak'))
        
        # test data
        self.assertListEqual(bak2.data.counts.tolist(), 
                             self.bak.data.counts.tolist())
        for i in range(bak2.num_chans):
            self.assertAlmostEqual(bak2.data.lo_edges[i],
                                   self.bak.data.lo_edges[i], places=4)
        for i in range(bak2.num_chans):
            self.assertAlmostEqual(bak2.data.hi_edges[i],
                                   self.bak.data.hi_edges[i], places=4)
        # test ebounds
        for i in range(bak2.num_chans):
            self.assertAlmostEqual(bak2.ebounds.low_edges()[i],
                                   self.bak.ebounds.low_edges()[i], places=4)
        for i in range(bak2.num_chans):
            self.assertAlmostEqual(bak2.ebounds.high_edges()[i],
                                   self.bak.ebounds.high_edges()[i], places=4)

        # test gti
        self.assertListEqual(bak2.gti.low_edges(), self.bak.gti.low_edges())
        #self.assertListEqual(bak2.gti.high_edges(), self.bak.gti.high_edges())

        # test attributes
        self.assertEqual(bak2.filename, 'test_bak.bak')
        self.assertEqual(bak2.num_chans, self.bak.num_chans)
        self.assertEqual(bak2.trigtime, self.bak.trigtime)
        
        # test headers
        self.assertEqual(bak2.headers['PRIMARY']['TSTART'], -899.0864419937134)
        self.assertEqual(bak2.headers['PRIMARY']['TSTOP'], -897.0384419937134)
        self.assertEqual(bak2.headers['PRIMARY']['TRIGTIME'], 356223561.133346)
        self.assertEqual(bak2.headers['EBOUNDS']['DETCHANS'], 8)
        self.assertEqual(bak2.headers['SPECTRUM']['DETCHANS'], 8)
        self.assertEqual(bak2.headers['SPECTRUM']['EXPOSURE'], 2.021)
        
        bak2.write(this_dir, overwrite=True)
        bak2.close()


if __name__ == '__main__':
    unittest.main()


