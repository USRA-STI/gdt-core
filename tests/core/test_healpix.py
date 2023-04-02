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
import os
import unittest
import healpy as hp

from gdt.core.healpix import *


class TestHealPix(unittest.TestCase):
    def setUp(self):
        hpx_arr = np.zeros(12)
        hpx_arr[0] = 0.5
        self.hpx = HealPix.from_data(hpx_arr, trigtime=1.0, filename='test.fit')

    def test_npix(self):
        self.assertEqual(self.hpx.npix, 12)

    def test_nside(self):
        self.assertEqual(self.hpx.nside, 1)
    
    def test_pixel_area(self):
        self.assertAlmostEqual(self.hpx.pixel_area, 3437.75, places=2)
    
    def test_trigtime(self):
        self.assertEqual(self.hpx.trigtime, 1.0)

    def test_multiply(self):
        hpx_mult = HealPix.multiply(self.hpx, self.hpx, output_nside=1)
        self.assertEqual(hpx_mult._hpx[0], 0.25)
        self.assertCountEqual(hpx_mult._hpx[1:], np.zeros(11))
        
        with self.assertRaises(TypeError):
            HealPix.multiply(self.hpx, 1.0)
        
        with self.assertRaises(TypeError):
            HealPix.multiply(1.0, self.hpx)

        with self.assertRaises(ValueError):
            HealPix.multiply(self.hpx, self.hpx, primary=2)

    def test_from_data_errors(self):
        with self.assertRaises(TypeError):
            HealPix.from_data(1.0)

        with self.assertRaises(TypeError):
            HealPix.from_data(self.hpx._hpx, trigtime='')

        with self.assertRaises(ValueError):
            HealPix.from_data(self.hpx._hpx, trigtime=-1.0)


class TestHealPixEffectiveArea(unittest.TestCase):
    def setUp(self):
        self.hpx = HealPixEffectiveArea.from_uniform(100.0, nside=2)
    
    def test_eff_area(self):
        self.assertCountEqual(self.hpx.eff_area, [100.0]*self.hpx.npix)        
    
    def test_effective_area(self):
        self.assertEqual(self.hpx.effective_area(0.0, 0.0), 100.0)
    
    def test_from_cosine(self):
        hpx = HealPixEffectiveArea.from_cosine(0.0, 0.0, 100.0)
        self.assertAlmostEqual(hpx.effective_area(0.0, 0.0), 100.0, delta=0.1)
        self.assertAlmostEqual(hpx.effective_area(0.0, 45.0), 70.7, delta=0.1)
        self.assertAlmostEqual(hpx.effective_area(0.0, 90.0), 0.0, delta=0.1)

        with self.assertRaises(TypeError):
            HealPixEffectiveArea.from_cosine('', 0.0, 100.0)
        with self.assertRaises(TypeError):
            HealPixEffectiveArea.from_cosine(0.0, '', 100.0)
        with self.assertRaises(TypeError):
            HealPixEffectiveArea.from_cosine(0.0, 0.0, '')
        with self.assertRaises(TypeError):
            HealPixEffectiveArea.from_cosine(0.0, 0.0, 100.0, coeff='')
        
        with self.assertRaises(ValueError):
            HealPixEffectiveArea.from_cosine(0.0, -5.0, 100.0)

        with self.assertRaises(ValueError):
            HealPixEffectiveArea.from_cosine(0.0, 0.0, -100.0)

        with self.assertRaises(ValueError):
            HealPixEffectiveArea.from_cosine(0.0, 0.0, 100.0, coeff=-1.0)
    
    def test_from_uniform_errors(self):
        with self.assertRaises(TypeError):
            HealPixEffectiveArea.from_uniform('')

        with self.assertRaises(ValueError):
            HealPixEffectiveArea.from_uniform(-1.0)
    
    def test_sum(self):
        sum_hpx = HealPixEffectiveArea.sum(self.hpx, self.hpx, output_nside=64)
        self.assertCountEqual(sum_hpx.eff_area, [200.0]*sum_hpx.npix)        
        
        with self.assertRaises(TypeError):
            HealPixEffectiveArea.sum(self.hpx, 1.0)
        
        with self.assertRaises(TypeError):
            HealPixEffectiveArea.sum(1.0, self.hpx)

        with self.assertRaises(ValueError):
            HealPixEffectiveArea.sum(self.hpx, self.hpx, primary=2)
    

class TestHealPixLocalization(unittest.TestCase):
    def setUp(self):
        self.hpx = HealPixLocalization.from_gaussian(0.0, 10.0, 10.0, nside=64)
    
    def test_centroid(self):
        self.assertAlmostEqual(self.hpx.centroid[0], 0.0, delta=0.2)
        self.assertAlmostEqual(self.hpx.centroid[1], 10.0, delta=0.2)
    
    def test_prob(self):
        self.assertEqual(len(self.hpx.prob), hp.nside2npix(64))
        self.assertGreaterEqual(self.hpx.prob.min(), 0.0)
        self.assertLessEqual(self.hpx.prob.max(), 1.0)

    def test_sig(self):
        self.assertEqual(len(self.hpx.sig), hp.nside2npix(64))
        self.assertGreaterEqual(self.hpx.sig.min(), 0.0)
        self.assertLessEqual(self.hpx.sig.max(), 1.0)

    def test_area(self):
        self.assertAlmostEqual(self.hpx.area(0.683), 715.0, delta=0.1)
        with self.assertRaises(ValueError):
            self.hpx.area(-0.1)

    def test_confidence(self):
        self.assertAlmostEqual(self.hpx.confidence(0.0, 10.0), 0.003, places=3)
        self.assertAlmostEqual(self.hpx.confidence(180.0, 10.0), 1.0, places=3)

    def test_confidence_region_path(self):
        path = self.hpx.confidence_region_path(0.01)
        for seg in path:
            diff = 360.0 % seg[(seg[:,0] > 0.0),0]
            for d in diff:
                self.assertLessEqual(d, 2.0)
            
            diff = np.abs(10.0-seg[:,1])
            for d in diff:
                self.assertLessEqual(d, 2.0)
        
        with self.assertRaises(ValueError):
            self.hpx.confidence_region_path(10.0)
        
    def test_convolve(self):
        def gauss_test():
            return (np.deg2rad([17.4]), [1.0])
        
        hpx = self.hpx.convolve(gauss_test)
        # this convolution results in ~2x increase in radius 
        # -> ~4x increase in area
        area_increase = hpx.area(0.5) / self.hpx.area(0.5)
        self.assertAlmostEqual(area_increase, 4.0, delta=0.1)

    def test_probability(self):
        self.assertAlmostEqual(self.hpx.probability(0.0, 10.0), 0.002, 
                               places=3)
        self.assertAlmostEqual(self.hpx.probability(180.0, 10.0), 0.0, 
                               places=9)

    def test_prob_array(self):
        prob_arr, ra_arr, dec_arr = self.hpx.prob_array()
        self.assertTupleEqual(prob_arr.shape, (180, 360))
        self.assertTupleEqual(ra_arr.shape, (360,))
        self.assertTupleEqual(dec_arr.shape, (180,))
        self.assertGreaterEqual(prob_arr.min(), 0.0)
        self.assertLessEqual(prob_arr.max(), 1.0)
        
        prob_arr, ra_arr, dec_arr = self.hpx.prob_array(sqdegrees=False, sig=True)
        self.assertGreaterEqual(prob_arr.min(), 0.0)
        self.assertLessEqual(prob_arr.max(), 1.0)

    def test_region_probability(self):
        prob = self.hpx.region_probability(self.hpx)
        self.assertGreaterEqual(prob, 0.97)
        
        hpx = HealPixLocalization.from_gaussian(180.0, 0.0, 1.0)
        prob = self.hpx.region_probability(hpx)
        self.assertLessEqual(prob, 0.001)
        
        with self.assertRaises(ValueError):
            self.hpx.region_probability(self.hpx, prior=-1.0)

    def test_source_probability(self):
        prob = self.hpx.source_probability(0.0, 10.0)
        self.assertGreaterEqual(prob, 0.97)

        prob = self.hpx.source_probability(180.0, 10.0)
        self.assertLessEqual(prob, 1e-6)

        with self.assertRaises(ValueError):
            self.hpx.region_probability(self.hpx, prior=-1.0)

    def test_from_annulus(self):
        hpx = HealPixLocalization.from_annulus(0.0, 0.0, 50.0, 5.0)
        
        with self.assertRaises(TypeError):
            HealPixLocalization.from_annulus('', 0.0, 50.0, 5.0)
        with self.assertRaises(TypeError):
            HealPixLocalization.from_annulus(0.0, '', 50.0, 5.0)
        with self.assertRaises(TypeError):
            HealPixLocalization.from_annulus(0.0, 0.0, '', 5.0)
        with self.assertRaises(TypeError):
            HealPixLocalization.from_annulus(0.0, 0.0, 50.0, '')
        
        with self.assertRaises(ValueError):
            HealPixLocalization.from_annulus(0.0, 100.0, 50.0, 5.0)
        with self.assertRaises(ValueError):
            HealPixLocalization.from_annulus(0.0, 0.0, -50.0, 5.0)
        with self.assertRaises(ValueError):
            HealPixLocalization.from_annulus(0.0, 0.0, 50.0, -5.0)

    def test_from_data(self):
        hpx = HealPixLocalization.from_data(self.hpx.prob) 

    def test_from_gaussian(self):
        with self.assertRaises(TypeError):
            HealPixLocalization.from_gaussian('', 0.0, 5.0)
        with self.assertRaises(TypeError):
            HealPixLocalization.from_gaussian(0.0, '', 5.0)
        with self.assertRaises(TypeError):
            HealPixLocalization.from_gaussian(0.0, 0.0, '')

        with self.assertRaises(ValueError):
            HealPixLocalization.from_gaussian(0.0, 100.0, 5.0)

        with self.assertRaises(ValueError):
            HealPixLocalization.from_gaussian(0.0, 0.0, -5.0)

    def test_from_vertices(self):
        ra_pts = [270.0, 180.0, 90.0]
        dec_pts = [15.0, 60.0, 15.0]
        hpx = HealPixLocalization.from_vertices(ra_pts, dec_pts)
        
        with self.assertRaises(ValueError):
            HealPixLocalization.from_vertices(ra_pts, [150.0, 60.0, 15.0])
        

if __name__ == '__main__':
    unittest.main()
      
