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
import unittest, os, shutil, sys
import numpy as np
from gdt.core.binning.binned import *


counts = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
exposure = np.array([1.024, 1.01, 1.00, 0.99, 1.02, 1.024, 0.80, 1.01])
edges = np.linspace(0.0, 8.192, 9)
index = np.array([0,2,4,6,8])

times = np.array([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 1.0, 1.01, 1.02, 1.03, 1.04])
    
class TestRebinningBinned(unittest.TestCase):
    
    def test_combine_by_factor(self): 
        new_counts, new_exposure, new_edges = combine_by_factor(counts, exposure, edges, 4)
        self.assertCountEqual(new_counts, np.array([10.0, 26.0]))
        self.assertCountEqual(new_exposure, np.array([4.024, 3.854]))
        self.assertCountEqual(new_edges, np.array([0.0, 4.096, 8.192]))
        
        self.assertRaises(AssertionError, combine_by_factor, counts, exposure, edges, 0)

    def test_rebin_by_time(self): 
        new_counts, new_exposure, new_edges = rebin_by_time(counts, exposure, edges, 4.096)
        self.assertCountEqual(new_counts, np.array([10.0, 26.0]))
        self.assertCountEqual(new_exposure, np.array([4.024, 3.854]))
        self.assertCountEqual(new_edges, np.array([0.0, 4.096, 8.192]))
        
        self.assertRaises(AssertionError, rebin_by_time, counts, exposure, edges, -1.0)

    def test_combine_into_one(self): 
        new_counts, new_exposure, new_edges = combine_into_one(counts, exposure, edges)
        self.assertCountEqual(new_counts, np.array([36.0]))
        self.assertCountEqual(new_exposure, np.array([7.878]))
        self.assertCountEqual(new_edges, np.array([0.0, 8.192]))
    
    def test_rebin_by_edge_index(self):
        new_counts, new_exposure, new_edges =rebin_by_edge_index(counts, exposure, edges, index)
        self.assertCountEqual(new_counts, np.array([3.0, 7.0, 11.0, 15.0]))
        self.assertCountEqual(new_exposure, np.array([2.034, 1.99, 2.044, 1.81]))
        self.assertCountEqual(new_edges, np.array([0.0, 2.048, 4.096, 6.144, 8.192]))
        
    def test_rebin_by_edge(self):
        new_edges = np.array([0.0, 2.048, 4.096, 6.144, 8.192])
        new_counts, new_exposure, new_edges =rebin_by_edge(counts, exposure, edges, new_edges)
        self.assertCountEqual(new_counts, np.array([3.0, 7.0, 11.0, 15.0]))
        self.assertCountEqual(new_exposure, np.array([2.034, 1.99, 2.044, 1.81]))
        self.assertCountEqual(new_edges, np.array([0.0, 2.048, 4.096, 6.144, 8.192]))

    def test_rebin_by_snr(self):
        back_counts = np.ones_like(counts)
        new_counts, new_exposure, new_edges = rebin_by_snr(counts, exposure, edges,
                                                           back_counts, 3.0)
        self.assertEqual(new_counts.size, 5)
        self.assertCountEqual(new_counts, np.array([10.0, 5.0, 6.0, 7.0, 8.0]))
        self.assertCountEqual(new_exposure, np.array([4.024, 1.02, 1.024, 0.8, 1.01]))
        self.assertCountEqual(new_edges, np.array([0.0, 4.096, 5.120, 6.144, 7.168, 8.192]))

        self.assertRaises(ValueError, rebin_by_snr, counts, exposure, edges, 
                          -1.0*back_counts, 3.0)


if __name__ == '__main__':
    unittest.main()
      
