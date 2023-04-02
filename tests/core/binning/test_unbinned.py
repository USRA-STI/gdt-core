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
from gdt.core.binning.unbinned import *


times = np.array([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 1.0, 1.01, 1.02, 1.03, 1.04])
    
class TestBinningUnbinned(unittest.TestCase):
    
    def test_bin_by_time(self):        
        new_edges = bin_by_time(times, 1.0)
        self.assertCountEqual(new_edges, np.array([0.0, 1.0, 2.0]))
        
        new_edges = bin_by_time(times, 0.5, time_ref=1.0)
        self.assertCountEqual(new_edges, np.array([0.0, 0.5, 1.0, 1.5]))
        
        self.assertRaises(AssertionError, bin_by_time, times, -1.0)
           
    def test_combine_into_one(self):
        new_edges = combine_into_one(times, 0.5, 2.0)
        self.assertCountEqual(new_edges, np.array([0.5, 2.0]))
        
        self.assertRaises(AssertionError, combine_into_one, times, 2.0, 0.5)
        self.assertRaises(ValueError, combine_into_one, times, 10.0, 20.0)
       
    def test_combine_by_factor(self): 
        old_edges = np.array([0.0, 0.5, 1.0, 1.5])
        new_edges = combine_by_factor(times, old_edges, 2)
        self.assertCountEqual(new_edges, np.array([0.0, 1.0]))

        self.assertRaises(AssertionError, combine_by_factor, times, old_edges, -1)
        
    def test_bin_by_snr(self):
        back_rates = np.ones_like(times)
        new_edges = bin_by_snr(times, back_rates, 5.0)
        
        self.assertCountEqual(new_edges, np.array([0.0, 0.02, 1.0, 1.01, 1.02, 1.03, 1.04]))
    
    def test_time_to_spill(self):
        
        new_edges = time_to_spill(times, 5)
        self.assertCountEqual(new_edges, np.array([0.0, 0.05, 1.04]))
        
        self.assertRaises(AssertionError, time_to_spill, times, -1)

    def test_bin_by_edges(self):
        edges = np.array([0.0, 0.25, 0.50, 0.75, 1.0])
        new_edges = bin_by_edges(times, edges)
        self.assertCountEqual(new_edges, edges)

if __name__ == '__main__':
    unittest.main()
      
