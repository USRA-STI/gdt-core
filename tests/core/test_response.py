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
from gdt.core.data_primitives import ResponseMatrix, Ebounds, EnergyBins, TimeBins
from gdt.core.headers import FileHeaders
from gdt.core.response import Rsp, Rsp2

class TestRsp(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        matrix = np.diag([0.1, 0.2, 0.3, 0.2, 0.1, 0.1])
        emin = [10., 20., 40., 80., 160., 320.]
        emax = [20., 40., 80., 160., 320., 640.]
        chanlo = [10.1, 20.1, 40.1, 80.1, 160.1, 320.1]
        chanhi = [20.1, 40.1, 80.1, 160.1, 320.1, 640.1]
        drm = ResponseMatrix(matrix, emin, emax, chanlo, chanhi)
        tstart = 356223432.1073
        tstop = 356223481.2601
        trigtime = 356223561.133346
        detector = 'det0'
        cls.rsp = Rsp.from_data(drm, start_time=tstart, stop_time=tstop,
                                 trigger_time=trigtime, detector=detector)

    def test_detector(self):
        self.assertEqual(self.rsp.detector, 'det0')
    
    def test_drm(self):
        self.assertIsInstance(self.rsp.drm, ResponseMatrix)

    def test_ebounds(self):
        self.assertIsInstance(self.rsp.ebounds, Ebounds)
        self.assertTupleEqual(self.rsp.ebounds.range, (10.1, 640.1))     
    
    def test_filename(self):
        self.assertIsNone(self.rsp.filename)    
    
    def test_headers(self):
        self.assertIsNone(self.rsp.headers)
    
    def test_num_chans(self):
        self.assertEqual(self.rsp.num_chans, 6)
    
    def test_num_ebins(self):
        self.assertEqual(self.rsp.num_ebins, 6)

    def test_tcent(self):
        self.assertAlmostEqual(self.rsp.tcent, -104.4496, places=4)
    
    def test_trigtime(self):
        self.assertEqual(self.rsp.trigtime, 356223561.133346)
    
    def test_tstart(self):
        self.assertAlmostEqual(self.rsp.tstart, -129.0260, places=4)

    def test_tstop(self):
        self.assertAlmostEqual(self.rsp.tstop, -79.8732, places=4)

    def test_fold_spectrum(self):
        def zeroth_poly(c, x):
            return np.full(x.shape, c[0])
        
        count_spec = self.rsp.fold_spectrum(zeroth_poly, (1.0,))
        self.assertIsInstance(count_spec, EnergyBins)
        self.assertListEqual(count_spec.counts.tolist(), 
                             [1., 4., 12., 16., 16., 32.])
        self.assertListEqual(count_spec.lo_edges.tolist(), 
                             [10.1, 20.1, 40.1, 80.1, 160.1, 320.1])
        self.assertListEqual(count_spec.hi_edges.tolist(), 
                             [20.1, 40.1, 80.1, 160.1, 320.1, 640.1])
        
        # non-unit exposure
        count_spec = self.rsp.fold_spectrum(zeroth_poly, (1.0,), exposure=2.0)
        self.assertListEqual(count_spec.counts.tolist(), 
                             [2., 8., 24., 32., 32., 64.])
        
        # channel mask
        mask = np.array([False, True, True, True, True, False])
        count_spec = self.rsp.fold_spectrum(zeroth_poly, (1.0,), 
                                            channel_mask=mask)
        self.assertEqual(count_spec.size, 4)
        self.assertListEqual(count_spec.counts.tolist(), [4., 12., 16., 16.])
        
        
        # invalid exposure type
        with self.assertRaises(TypeError):
            self.rsp.fold_spectrum(zeroth_poly, (1.0,), exposure='')

        # invalid exposure value
        with self.assertRaises(ValueError):
            self.rsp.fold_spectrum(zeroth_poly, (1.0,), exposure=-1.0)

    def test_rebin(self):
        rebinned = self.rsp.rebin(factor=2)
        self.assertEqual(rebinned.num_chans, 3)
        self.assertEqual(rebinned.num_ebins, 6)
        self.assertListEqual(rebinned.ebounds.low_edges(), [10.1, 40.1, 160.1])
        self.assertListEqual(rebinned.ebounds.high_edges(), [40.1, 160.1, 640.1])
        self.assertAlmostEqual(rebinned.tstart, -129.0260, places=4)
        self.assertAlmostEqual(rebinned.tstop, -79.8732, places=4)
        self.assertEqual(rebinned.trigtime, 356223561.133346)

    def test_resample(self):
        resampled = self.rsp.resample(num_photon_bins=3)
        self.assertEqual(resampled.num_chans, 6)
        self.assertEqual(resampled.num_ebins, 3)
        self.assertAlmostEqual(resampled.tstart, -129.0260, places=4)
        self.assertAlmostEqual(resampled.tstop, -79.8732, places=4)
        self.assertEqual(resampled.trigtime, 356223561.133346)

    def test_write(self):
        with self.assertRaises(NameError):
            self.rsp.write('.')

    def test_errors(self):
        with self.assertRaises(TypeError):
            Rsp.from_data(self.rsp.ebounds)
        
        with self.assertRaises(ValueError):
            Rsp.from_data(self.rsp.drm, trigger_time=-1.0)
        
        with self.assertRaises(TypeError):
            Rsp.from_data(self.rsp.drm, headers=1.0)
 

class TestRsp2(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        emin = [10., 20., 40., 80., 160., 320.]
        emax = [20., 40., 80., 160., 320., 640.]
        chanlo = [10.1, 20.1, 40.1, 80.1, 160.1, 320.1]
        chanhi = [20.1, 40.1, 80.1, 160.1, 320.1, 640.1]
        trigtime = 356223561.133346
        detector = 'det0'

        matrix = np.diag([0.1, 0.2, 0.3, 0.2, 0.1, 0.1])
        drm1 = ResponseMatrix(matrix, emin, emax, chanlo, chanhi)
        tstart1 = 356223432.1073
        tstop1 = 356223481.2601
        rsp1 = Rsp.from_data(drm1, start_time=tstart1, stop_time=tstop1,
                             trigger_time=trigtime, detector=detector)

        matrix = np.diag([0.1, 0.2, 0.3, 0.2, 0.1, 0.0])
        drm2 = ResponseMatrix(matrix, emin, emax, chanlo, chanhi)
        tstart2 = 356223481.2601
        tstop2 = 356223530.4129
        rsp2 = Rsp.from_data(drm2, start_time=tstart2, stop_time=tstop2,
                             trigger_time=trigtime, detector=detector)

        matrix = np.diag([0.1, 0.2, 0.3, 0.2, 0.0, 0.0])
        drm3 = ResponseMatrix(matrix, emin, emax, chanlo, chanhi)
        tstart3 = 356223530.4129
        tstop3 = 356223575.4696
        rsp3 = Rsp.from_data(drm3, start_time=tstart3, stop_time=tstop3,
                             trigger_time=trigtime, detector=detector)
        
        cls.rsp2 = Rsp2.from_rsps([rsp1, rsp2, rsp3])

    def test_detector(self):
        self.assertEqual(self.rsp2.detector, 'det0')

    def test_ebounds(self):
        self.assertIsInstance(self.rsp2.ebounds, Ebounds)
        self.assertTupleEqual(self.rsp2.ebounds.range, (10.1, 640.1))     

    def test_filename(self):
        self.assertIsNone(self.rsp2.filename)    
    
    def test_num_chans(self):
        self.assertEqual(self.rsp2.num_chans, 6)

    def test_num_drms(self):
        self.assertEqual(self.rsp2.num_drms, 3)
    
    def test_num_ebins(self):
        self.assertEqual(self.rsp2.num_ebins, 6)

    def test_tcent(self):
        tcent = [-104.4496, -55.2968, -8.1921]
        for i in range(self.rsp2.num_drms):
            self.assertAlmostEqual(self.rsp2.tcent[i], tcent[i], places=4)
 
    def test_trigtime(self):
        self.assertEqual(self.rsp2.trigtime, 356223561.133346)
    
    def test_tstart(self):
        tstart = [-129.0260, -79.8732, -30.7204]
        for i in range(self.rsp2.num_drms):
            self.assertAlmostEqual(self.rsp2.tstart[i], tstart[i], places=4)

    def test_tstop(self):
        tstop = [-79.8732, -30.7204, 14.3363]
        for i in range(self.rsp2.num_drms):
            self.assertAlmostEqual(self.rsp2.tstop[i], tstop[i], places=4)

    def test_drm_index(self):
        # one index
        self.assertListEqual(self.rsp2.drm_index((-100.0, -80.0)).tolist(), [0])

        # multiple indices
        self.assertListEqual(self.rsp2.drm_index((-100.0, 0.0)).tolist(), 
                             [0, 1, 2])
        
        # bad range
        with self.assertRaises(AssertionError):
            self.rsp2.drm_index((0.0, -10.0))
    
    def test_extract_drm(self):
        # extract by index
        rsp = self.rsp2.extract_drm(index=1)
        self.assertAlmostEqual(rsp.tstart, -79.8732, places=4)
        self.assertAlmostEqual(rsp.tstop, -30.7204, places=4)
        
        # extract by time
        rsp = self.rsp2.extract_drm(atime=-50.0)
        self.assertAlmostEqual(rsp.tstart, -79.8732, places=4)
        self.assertAlmostEqual(rsp.tstop, -30.7204, places=4)
        
        # neither index nor atime is set
        with self.assertRaises(RuntimeError):
            self.rsp2.extract_drm()        

    def test_interpolate(self):
        # at the tcent of a DRM
        interp = self.rsp2.interpolate(-104.4496)
        effarea = [0.1, 0.2, 0.3, 0.2, 0.1, 0.1]
        for i in range(6):
            self.assertAlmostEqual(interp.drm.matrix[i,i], effarea[i], places=1)
        self.assertAlmostEqual(interp.tstart, -104.4496, places=4)
        self.assertAlmostEqual(interp.tstop, -104.4496, places=4)
        self.assertEqual(interp.trigtime, 356223561.133346)
        
        # halfway between tcents of two DRMs
        interp = self.rsp2.interpolate(-79.8732)
        effarea = [0.1, 0.2, 0.3, 0.2, 0.1, 0.05]
        for i in range(6):
            self.assertAlmostEqual(interp.drm.matrix[i,i], effarea[i], places=2)
        
        # only a single DRM RSP2
        rsp2 = Rsp2.from_rsps([self.rsp2[0]])
        with self.assertRaises(ValueError):
            rsp2.interpolate(-79.8732)

    def test_nearest_drm(self):
        # before the first DRM
        rsp = self.rsp2.nearest_drm(-200.0)
        self.assertAlmostEqual(rsp.tstart, -129.0260, places=4)
        self.assertAlmostEqual(rsp.tstop, -79.8732, places=4)

        # in the DRM range
        rsp = self.rsp2.nearest_drm(-50.0)
        self.assertAlmostEqual(rsp.tstart, -79.8732, places=4)
        self.assertAlmostEqual(rsp.tstop, -30.7204, places=4)

        # after the last DRM
        rsp = self.rsp2.nearest_drm(50.0)
        self.assertAlmostEqual(rsp.tstart, -30.7204, places=4)
        self.assertAlmostEqual(rsp.tstop, 14.3363, places=4)
    
    def test_weighted(self):
        counts = [1, 2, 3, 4]
        tstart = [-129.0260, -104.4496, -79.8732, -55.2968]
        tstop = [-104.4496, -79.8732, -55.2968, -30.7204]
        exposure = [24.576] * 4
        bins = TimeBins(counts, tstart, tstop, exposure)
        
        # nearest DRM
        weighted = self.rsp2.weighted(bins)
        effarea = [0.1, 0.2, 0.3, 0.2, 0.1, 0.03]
        for i in range(6):
            self.assertAlmostEqual(weighted.drm.matrix[i,i], effarea[i], places=2)
        self.assertAlmostEqual(weighted.tstart, -129.0260, places=4)
        self.assertAlmostEqual(weighted.tstop, -30.7204, places=4)
        self.assertEqual(weighted.trigtime, 356223561.133346)

        # interpolated        
        weighted = self.rsp2.weighted(bins, interpolate=True)
        effarea = [0.1, 0.2, 0.3, 0.2, 0.0896, 0.0325]
        for i in range(6):
            self.assertAlmostEqual(weighted.drm.matrix[i,i], effarea[i], places=4)

        # only a single DRM RSP2
        rsp2 = Rsp2.from_rsps([self.rsp2[0]])
        weighted = rsp2.weighted(bins)
        effarea = [0.1, 0.2, 0.3, 0.2, 0.1, 0.1]
        for i in range(6):
            self.assertAlmostEqual(weighted.drm.matrix[i,i], effarea[i], places=2)
        
        # incorrect input
        with self.assertRaises(TypeError):
            self.rsp2.weighted(1.0)

    def test_write(self):
        with self.assertRaises(NameError):
            self.rsp2.write('.')

    def test_index_access(self):
        rsp = self.rsp2[1]
        self.assertAlmostEqual(rsp.tstart, -79.8732, places=4)
        self.assertAlmostEqual(rsp.tstop, -30.7204, places=4)
        
        # invalid index
        with self.assertRaises(IndexError):
            self.rsp2[4]

    def test_errors(self):
        # not a list of Rsp
        with self.assertRaises(TypeError):
            Rsp2.from_rsps(self.rsp2[0])
        
        # one in the list is not an Rsp
        with self.assertRaises(TypeError):
            Rsp2.from_rsps([self.rsp2[0], 1.0])
                           
if __name__ == '__main__':
    unittest.main()


