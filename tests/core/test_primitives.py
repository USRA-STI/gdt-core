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
import unittest
import numpy as np
from gdt.core.data_primitives import *
from gdt.core.binning.binned import combine_by_factor
from gdt.core.binning.unbinned import bin_by_time

class TestRange(unittest.TestCase):
    
    def setUp(self):
        self.range = Range(0.0, 10.0)
    
    def test_center(self):
        self.assertEqual(self.range.center, 5.0)
    
    def test_width(self):
        self.assertEqual(self.range.width, 10.0)
    
    def test_as_tuple(self):
        self.assertTupleEqual(self.range.as_tuple(), (0.0, 10.0))
    
    def test_contains(self):
        self.assertTrue(self.range.contains(1.0, inclusive=False))
        self.assertFalse(self.range.contains(10.0, inclusive=False))
        self.assertFalse(self.range.contains(0.0, inclusive=False))
        self.assertFalse(self.range.contains(20.0, inclusive=False))

        self.assertTrue(self.range.contains(1.0, inclusive=True))
        self.assertTrue(self.range.contains(10.0, inclusive=True))
        self.assertTrue(self.range.contains(0.0, inclusive=True))
        self.assertFalse(self.range.contains(20.0, inclusive=False))
    
    def test_intersection(self):
        # intersect on high end
        range2 = Range(5.0, 15.0)
        range_intersect = Range.intersection(self.range, range2)
        self.assertTupleEqual(range_intersect.as_tuple(), (5.0, 10.0))
        
        # flip order
        range_intersect = Range.intersection(range2, self.range)
        self.assertTupleEqual(range_intersect.as_tuple(), (5.0, 10.0))
        
        # intersect on low end
        range2 = Range(-5.0, 5.0)
        range_intersect = Range.intersection(self.range, range2)
        self.assertTupleEqual(range_intersect.as_tuple(), (0.0, 5.0))
        
        # one range contained inside other
        range2 = Range(5.0, 6.0)
        range_intersect = Range.intersection(self.range, range2)
        self.assertTupleEqual(range_intersect.as_tuple(), (5.0, 6.0))
        
        # no intersection
        range2 = Range(20.0, 30.0)
        range_intersect = Range.intersection(self.range, range2)
        self.assertIsNone(range_intersect)
        
    def test_union(self):
        # union on high end
        range2 = Range(5.0, 15.0)
        range_union = Range.union(self.range, range2)        
        self.assertTupleEqual(range_union.as_tuple(), (0.0, 15.0))
        
        # flip order
        range_union = Range.union(range2, self.range)
        self.assertTupleEqual(range_union.as_tuple(), (0.0, 15.0))

        # union on low end
        range2 = Range(-5.0, 5.0)
        range_union = Range.union(self.range, range2)
        self.assertTupleEqual(range_union.as_tuple(), (-5.0, 10.0))

        # one range contained inside other
        range2 = Range(5.0, 6.0)
        range_union = Range.union(self.range, range2)
        self.assertTupleEqual(range_union.as_tuple(), (0.0, 10.0))

    def test_equal(self):
        self.assertTrue(self.range == self.range)
        
        range2 = Range(1.0, 10.0)
        self.assertFalse(range2 == self.range)
        self.assertFalse(self.range == range2)

    def test_translate_negative(self):
        range2 = self.range.translate(-5.0)
        self.assertTupleEqual(range2.as_tuple(), (-5.0, 5.0))

    def test_translate_positive(self):
        range2 = self.range.translate(5.0)
        self.assertTupleEqual(range2.as_tuple(), (5.0, 15.0))
        

class TestTimeRange(unittest.TestCase):
    
    def setUp(self):
        self.range = TimeRange(-5.0, 100.0)
    
    def test_duration(self):
        self.assertEqual(self.range.duration, 105.0)
    
    def test_tstart(self):
        self.assertEqual(self.range.tstart, -5.0)

    def test_tstop(self):
        self.assertEqual(self.range.tstop, 100.0)
    
    def test_errors(self):
        with self.assertRaises(TypeError):
            range = TimeRange('', 1.0)
        with self.assertRaises(TypeError):
            range = TimeRange(1.0, '')


class TestEnergyRange(unittest.TestCase):
    
    def setUp(self):
        self.range = EnergyRange(10.0, 1000.0)
    
    def test_emax(self):
        self.assertEqual(self.range.emax, 1000.0)

    def test_emin(self):
        self.assertEqual(self.range.emin, 10.0)
    
    def test_log_center(self):
        self.assertEqual(self.range.log_center, 100.0)
    
    def test_errors(self):
        with self.assertRaises(TypeError):
            range = EnergyRange('', 1.0)
        with self.assertRaises(TypeError):
            range = EnergyRange(1.0, '')
        

class TestIntervals(unittest.TestCase):
    
    def setUp(self):
        low = (0.0, 10.0, 20.0)
        high = (10.0, 20.0, 30.0)
        self.intervals = Intervals.from_bounds(low, high)
    
    def test_init(self):
        intervals = Intervals(interval=Range(0.0, 10.0))
        self.assertTupleEqual(intervals[0].as_tuple(), (0.0, 10.0))
    
    def test_intervals(self):
        self.assertTupleEqual(self.intervals.intervals[0].as_tuple(), (0.0, 10.0))
        self.assertTupleEqual(self.intervals.intervals[1].as_tuple(), (10.0, 20.0))
        self.assertTupleEqual(self.intervals.intervals[2].as_tuple(), (20.0, 30.0))
    
    def test_access(self):
        self.assertTupleEqual(self.intervals[0].as_tuple(), (0.0, 10.0))
        self.assertTupleEqual(self.intervals[1].as_tuple(), (10.0, 20.0))
        self.assertTupleEqual(self.intervals[2].as_tuple(), (20.0, 30.0))
    
    def test_num_intervals(self):
        self.assertEqual(self.intervals.num_intervals, 3)
    
    def test_range(self):
        self.assertTupleEqual(self.intervals.range, (0.0, 30.0))
    
    def test_as_list(self):
        vals = [(0.0, 10.0), (10.0, 20.0), (20.0, 30.0)]
        lst = self.intervals.as_list()
        for i in range(3):
            self.assertTupleEqual(lst[i], vals[i])

    def test_contains(self):
        self.assertTrue(self.intervals.contains(1.0, inclusive=False))
        self.assertFalse(self.intervals.contains(10.0, inclusive=False))
        self.assertFalse(self.intervals.contains(35.0, inclusive=False))

        self.assertTrue(self.intervals.contains(1.0, inclusive=True))
        self.assertTrue(self.intervals.contains(10.0, inclusive=True))
        self.assertFalse(self.intervals.contains(35.0, inclusive=True))

    def test_high_edges(self):
        self.assertListEqual(self.intervals.high_edges(), [10.0, 20.0, 30.0])

    def test_index(self):
        assert self.intervals.index(5.0) == 0
        assert self.intervals.index(20.0) == 1
        assert self.intervals.index(-5.0) is None
        assert self.intervals.index(50.0) is None

    def test_insert(self):
    
        # insert a duplicate
        self.intervals.insert(Range(0.0, 10.0))
        vals = [(0.0, 10.0), (10.0, 20.0), (20.0, 30.0)]
        lst = self.intervals.as_list()
        for i in range(3):
            self.assertTupleEqual(lst[i], vals[i])

        # insert at end
        self.intervals.insert(Range(30.0, 40.0))
        vals = [(0.0, 10.0), (10.0, 20.0), (20.0, 30.0), (30.0, 40.0)]
        lst = self.intervals.as_list()
        for i in range(4):
            self.assertTupleEqual(lst[i], vals[i])
        
        # insert at beginning
        self.intervals.insert(Range(-10.0, -1.0))
        vals = [(-10.0, -1.0), (0.0, 10.0), (10.0, 20.0), (20.0, 30.0), 
                (30.0, 40.0)]
        lst = self.intervals.as_list()
        for i in range(5):
            self.assertTupleEqual(lst[i], vals[i])
        
        # insert in middle and overlapping
        self.intervals.insert(Range(15.0, 25.0))
        vals = [(-10.0, -1.0), (0.0, 10.0), (10.0, 30.0), (30.0, 40.0)]
        lst = self.intervals.as_list()
        for i in range(4):
            self.assertTupleEqual(lst[i], vals[i])
                
    def test_low_edges(self):
        self.assertListEqual(self.intervals.low_edges(), [0.0, 10.0, 20.0])
    
    def test_from_list(self):
        lst = [(0.0, 10.0), (10.0, 20.0), (20.0, 30.0)]
        intervals = Intervals.from_list(lst)
        self.assertListEqual(intervals.low_edges(), [0.0, 10.0, 20.0])
        self.assertListEqual(intervals.high_edges(), [10.0, 20.0, 30.0])
    
    def test_intersection(self):
        # overlapping
        intervals2 = Intervals.from_list([(-5.0, 5.0), (25.0, 35.0)])
        intervals_intersect = Intervals.intersection(self.intervals, intervals2)
        self.assertListEqual(intervals_intersect.high_edges(), [5.0, 30.0])
        self.assertListEqual(intervals_intersect.low_edges(), [0.0, 25.0])

        # flip order
        intervals_intersect = Intervals.intersection(intervals2, self.intervals)
        self.assertListEqual(intervals_intersect.high_edges(), [5.0, 30.0])
        self.assertListEqual(intervals_intersect.low_edges(), [0.0, 25.0])

        # one encompasses the other
        intervals2 = Intervals(Range(-10.0, 50.0))
        intervals_intersect = Intervals.intersection(self.intervals, intervals2)
        self.assertListEqual(intervals_intersect.high_edges(), [10.0, 20.0, 30.0])
        self.assertListEqual(intervals_intersect.low_edges(), [0.0, 10.0, 20.0])
        
        # no overlap
        intervals2 = Intervals(Range(50.0, 100.0))
        intervals_intersect = Intervals.intersection(self.intervals, intervals2)
        self.assertListEqual(intervals_intersect.high_edges(), [])
        self.assertListEqual(intervals_intersect.low_edges(), [])

    def test_merge(self):
        # overlapping
        intervals2 = Intervals.from_list([(-5.0, 5.0), (25.0, 35.0)])
        intervals_merged = Intervals.merge(self.intervals, intervals2)
        self.assertListEqual(intervals_merged.high_edges(), [10.0, 20.0, 35.0])
        self.assertListEqual(intervals_merged.low_edges(), [-5.0, 10.0, 20.0])

        # flip order
        intervals_merged = Intervals.merge(intervals2, self.intervals)
        self.assertListEqual(intervals_merged.high_edges(), [10.0, 20.0, 35.0])
        self.assertListEqual(intervals_merged.low_edges(), [-5.0, 10.0, 20.0])

        # one encompasses the other
        intervals2 = Intervals(Range(-10.0, 50.0))
        intervals_merged = Intervals.merge(self.intervals, intervals2)
        self.assertListEqual(intervals_merged.high_edges(), [50.0])
        self.assertListEqual(intervals_merged.low_edges(), [-10.0])

        # no overlap
        intervals2 = Intervals(Range(50.0, 100.0))
        intervals_merged = Intervals.merge(self.intervals, intervals2)
        self.assertListEqual(intervals_merged.high_edges(), [10.0, 20.0, 30.0, 100.0])
        self.assertListEqual(intervals_merged.low_edges(), [0.0, 10.0, 20.0, 50.0])


class TestGti(unittest.TestCase):
    
    def test_init(self):
        gti = Gti.from_list([(0.0, 10.0)])
        self.assertIsInstance(gti[0], TimeRange)
    
    def test_from_boolean_mask(self):
        times = np.arange(10.0)

        # all times valid
        mask = np.ones_like(times, dtype=bool)
        gti = Gti.from_boolean_mask(times, mask)
        self.assertListEqual(gti.high_edges(), [9.0])
        self.assertListEqual(gti.low_edges(), [0.0])
        
        # exclude first
        mask[0] = False
        gti = Gti.from_boolean_mask(times, mask)
        self.assertListEqual(gti.high_edges(), [9.0])
        self.assertListEqual(gti.low_edges(), [1.0])
        
        # exclude last
        mask[-1] = False
        gti = Gti.from_boolean_mask(times, mask)
        self.assertListEqual(gti.high_edges(), [8.0])
        self.assertListEqual(gti.low_edges(), [1.0])

        # exclude in middle
        mask[4:6] = False
        gti = Gti.from_boolean_mask(times, mask)
        self.assertListEqual(gti.high_edges(), [3.0, 8.0])
        self.assertListEqual(gti.low_edges(), [1.0, 6.0])
        
        # no times valid
        mask[:] = False
        gti = Gti.from_boolean_mask(times, mask)
        self.assertIsNone(gti)
    
        
class TestEbounds(unittest.TestCase):
    
    def test_init(self):
        ebounds = Ebounds.from_list([(50.0, 300.0)])
        self.assertIsInstance(ebounds[0], EnergyRange)
    

class TestBins(unittest.TestCase):
    
    def setUp(self):
        counts = [0, 9, 16, 4]
        low = [0.0, 1.0, 2.0, 3.0]
        high = [1.0, 2.0, 3.0, 4.0]
        self.bins = Bins(counts, low, high)
    
    def test_centroids(self):
        self.assertListEqual(self.bins.centroids.tolist(), [0.5, 1.5, 2.5, 3.5])

    def test_counts(self):
        self.assertListEqual(self.bins.counts.tolist(), [0, 9, 16, 4])

    def test_count_uncertainty(self):
        self.assertListEqual(self.bins.count_uncertainty.tolist(), [0, 3, 4, 2])

    def test_hi_edges(self):
        self.assertListEqual(self.bins.hi_edges.tolist(), [1.0, 2.0, 3.0, 4.0])

    def test_lo_edges(self):
        self.assertListEqual(self.bins.lo_edges.tolist(), [0.0, 1.0, 2.0, 3.0])

    def test_range(self):
        self.assertTupleEqual(self.bins.range, (0.0, 4.0))

    def test_rates(self):
        self.assertListEqual(self.bins.rates.tolist(), [0.0, 9.0, 16.0, 4.0])

    def test_rate_uncertainty(self):
        self.assertListEqual(self.bins.rate_uncertainty.tolist(), 
                             [0.0, 3.0, 4.0, 2.0])

    def test_size(self):
        self.assertEqual(self.bins.size, 4)
    
    def test_widths(self):
        self.assertListEqual(self.bins.widths.tolist(), [1.0, 1.0, 1.0, 1.0])

    def test_closest_edge(self):
        # closest low edge
        self.assertEqual(self.bins.closest_edge(0.7, which='low'), 0.0)
        # closest high edge
        self.assertEqual(self.bins.closest_edge(0.3, which='high'), 1.0)
        # closest edge
        self.assertEqual(self.bins.closest_edge(0.3, which='either'), 0.0)

    def test_slice(self):
        # middle slice       
        bins2 = self.bins.slice(1.3, 2.3)
        self.assertTupleEqual(bins2.range, (1.0, 3.0))
        
        # slice below lower boundary
        bins2 = self.bins.slice(-1.0, 2.3)
        self.assertTupleEqual(bins2.range, (0.0, 3.0))

        # slice above upper boundary
        bins2 = self.bins.slice(1.3, 5.0)
        self.assertTupleEqual(bins2.range, (1.0, 4.0))
        
        # slice covering full range
        bins2 = self.bins.slice(-1.0, 5.0)
        self.assertTupleEqual(bins2.range, (0.0, 4.0))

        # slice one bin
        bins2 = self.bins.slice(2.1, 2.1)
        self.assertTupleEqual(bins2.range, (2.0, 3.0))

        # slice fully outside range
        bins2 = self.bins.slice(-2.0, -1.0)
        self.assertIsNone(bins2.range)
        
    def test_init_errors(self):
        with self.assertRaises(TypeError):
            Bins(1.0, self.bins.lo_edges, self.bins.hi_edges)
        
        with self.assertRaises(TypeError):
            Bins(self.bins.counts, 1.0, self.bins.hi_edges)

        with self.assertRaises(TypeError):
            Bins(self.bins.counts, self.bins.lo_edges, 1.0)
        
        with self.assertRaises(ValueError):
            Bins(self.bins.counts[1:], self.bins.lo_edges, self.bins.hi_edges)


class TestExposureBins(unittest.TestCase):
    
    def setUp(self):
        counts = [0, 20, 12, 3, 20, 12]
        exposure = [0.9, 0.8, 0.75, 0.75, 0.8, 0.8]
        low = [0.0, 1.0, 2.5, 3.5, 4.5, 5.5]
        high = [1.0, 2.0, 3.5, 4.5, 5.5, 6.5]
        self.bins = ExposureBins(counts, low, high, exposure)

    def test_exposure(self):
        self.assertListEqual(self.bins.exposure.tolist(), 
                             [0.9, 0.8, 0.75, 0.75, 0.8, 0.8])

    def test_rates(self):
        self.assertListEqual(self.bins.rates.tolist(), 
                             [0.0, 25.0, 16.0, 4.0, 25.0, 15.0])

    def test_rate_uncertainty(self):
        uncert = self.bins.rate_uncertainty.tolist()
        test = [0.0, 5.59, 4.62, 2.31, 5.59, 4.33]
        for i in range(4):
            self.assertAlmostEqual(uncert[i], test[i], places=2)
    
    def test_contiguous_bins(self):
        cont_bins = self.bins.contiguous_bins()
        self.assertEqual(len(cont_bins), 2)
        self.assertTupleEqual(cont_bins[0].range, (0.0, 2.0))
        self.assertTupleEqual(cont_bins[1].range, (2.5, 6.5))
        
        self.assertEqual(len(cont_bins[0].contiguous_bins()), 1)
    
    def test_rebin(self):

        # rebin full range
        rebinned = self.bins.rebin(combine_by_factor, 2)
        self.assertEqual(rebinned.size, 3)
        self.assertListEqual(rebinned.hi_edges.tolist(), [2.0, 4.5, 6.5])
        self.assertListEqual(rebinned.lo_edges.tolist(), [0.0, 2.5, 4.5])
        self.assertListEqual(rebinned.counts.tolist(), [20, 15, 32])
        exp = [1.7, 1.5, 1.6]
        for i in range(3):
            self.assertAlmostEqual(rebinned.exposure[i], exp[i], places=2)

        # rebin from tstart through end of range        
        rebinned = self.bins.rebin(combine_by_factor, 2, tstart=2.7)
        self.assertEqual(rebinned.size, 4)
        self.assertListEqual(rebinned.hi_edges.tolist(), [1.0, 2.0, 4.5, 6.5])
        self.assertListEqual(rebinned.lo_edges.tolist(), [0.0, 1.0, 2.5, 4.5])
        self.assertListEqual(rebinned.counts.tolist(), [0, 20, 15, 32])
        exp = [0.9, 0.8, 1.5, 1.6]
        for i in range(4):
            self.assertAlmostEqual(rebinned.exposure[i], exp[i], places=2)

        # rebin from beginning of range through tstop       
        rebinned = self.bins.rebin(combine_by_factor, 2, tstop=1.7)
        self.assertEqual(rebinned.size, 5)
        self.assertListEqual(rebinned.hi_edges.tolist(), [2.0, 3.5, 4.5, 5.5, 6.5])
        self.assertListEqual(rebinned.lo_edges.tolist(), [0.0, 2.5, 3.5, 4.5, 5.5])
        self.assertListEqual(rebinned.counts.tolist(), [20, 12, 3, 20, 12])
        exp = [1.7, 0.75, 0.75, 0.8, 0.8]
        for i in range(5):
            self.assertAlmostEqual(rebinned.exposure[i], exp[i], places=2)

        # rebin middle      
        rebinned = self.bins.rebin(combine_by_factor, 2, tstart=0.5, tstop=3.7)
        self.assertEqual(rebinned.size, 4)
        self.assertListEqual(rebinned.hi_edges.tolist(), [2.0, 4.5, 5.5, 6.5])
        self.assertListEqual(rebinned.lo_edges.tolist(), [0.0, 2.5, 4.5, 5.5])
        self.assertListEqual(rebinned.counts.tolist(), [20, 15, 20, 12])
        exp = [1.7, 1.5, 0.8, 0.8]
        for i in range(4):
            self.assertAlmostEqual(rebinned.exposure[i], exp[i], places=2)

        # rebin outside range
        rebinned = self.bins.rebin(combine_by_factor, 2, tstart=7.0, tstop=10.)
        self.assertEqual(rebinned.size, 6)
        self.assertListEqual(rebinned.hi_edges.tolist(), [1.0, 2.0, 3.5, 4.5, 5.5, 6.5])
        self.assertListEqual(rebinned.lo_edges.tolist(), [0.0, 1.0, 2.5, 3.5, 4.5, 5.5])
        self.assertListEqual(rebinned.counts.tolist(), [0, 20, 12, 3, 20, 12])
        exp = [0.9, 0.8, 0.75, 0.75, 0.8, 0.8]
        for i in range(6):
            self.assertAlmostEqual(rebinned.exposure[i], exp[i], places=2)

    def test_slice(self):
        # middle slice       
        bins2 = self.bins.slice(1.3, 2.7)
        self.assertTupleEqual(bins2.range, (1.0, 3.5))
        
        # slice below lower boundary
        bins2 = self.bins.slice(-1.0, 2.7)
        self.assertTupleEqual(bins2.range, (0.0, 3.5))

        # slice above upper boundary
        bins2 = self.bins.slice(1.3, 7.0)
        self.assertTupleEqual(bins2.range, (1.0, 6.5))
        
        # slice covering full range
        bins2 = self.bins.slice(-1.0, 7.0)
        self.assertTupleEqual(bins2.range, (0.0, 6.5))

        # slice one bin
        bins2 = self.bins.slice(1.9, 1.9)
        self.assertTupleEqual(bins2.range, (1.0, 2.0))

        # slice fully outside range
        bins2 = self.bins.slice(-2.0, -1.0)
        self.assertIsNone(bins2.range)
    
    def test_merge(self):
        
        counts = [5, 7, 8]
        exposure = [0.9, 0.9, 0.9]
        
        # merge at high end
        low = [6.5, 7.5, 8.5]
        high = [7.5, 8.5, 9.5]
        bins2 = ExposureBins(counts, low, high, exposure)
        bins_merged = ExposureBins.merge([self.bins, bins2])
        self.assertListEqual(bins_merged.hi_edges.tolist(), [1.0, 2.0, 3.5, 4.5,
                                                             5.5, 6.5, 7.5, 8.5,
                                                             9.5])
        self.assertListEqual(bins_merged.lo_edges.tolist(), [0.0, 1.0, 2.5, 3.5,
                                                             4.5, 5.5, 6.5, 7.5,
                                                             8.5])

        
        # flip order
        bins_merged = ExposureBins.merge([bins2, self.bins])
        self.assertListEqual(bins_merged.hi_edges.tolist(), [1.0, 2.0, 3.5, 4.5,
                                                             5.5, 6.5, 7.5, 8.5,
                                                             9.5])
        self.assertListEqual(bins_merged.lo_edges.tolist(), [0.0, 1.0, 2.5, 3.5,
                                                             4.5, 5.5, 6.5, 7.5,
                                                             8.5])
        
        # merge at low end
        low = [-3.0, -2.0, -1.0]
        high = [-2.0, -1.0, 0.0]
        bins2 = ExposureBins(counts, low, high, exposure)
        bins_merged = ExposureBins.merge([self.bins, bins2])
        self.assertListEqual(bins_merged.hi_edges.tolist(), [-2.0, -1.0, 0.0, 
                                                             1.0, 2.0, 3.5, 4.5,
                                                             5.5, 6.5])
        self.assertListEqual(bins_merged.lo_edges.tolist(), [-3.0, -2.0, -1.0, 
                                                             0.0, 1.0, 2.5, 3.5,
                                                             4.5, 5.5])
        
        # overlapping merge (not allowed)
        low = [-1.5, -0.5, 0.5]
        high = [-0.5, 0.5, 1.5]
        bins2 = ExposureBins(counts, low, high, exposure)
        with self.assertRaises(ValueError):
            bins_merged = ExposureBins.merge([self.bins, bins2])
        with self.assertRaises(ValueError):
            bins_merged = ExposureBins.merge([bins2, self.bins])

    def test_sum(self):
        
        # same exposure
        bins2 = ExposureBins(self.bins.counts, self.bins.lo_edges, 
                             self.bins.hi_edges, self.bins.exposure)
        bins_summed = ExposureBins.sum([self.bins, bins2])
        self.assertListEqual(bins_summed.counts.tolist(), [0, 40, 24, 6, 40, 24])
        self.assertListEqual(bins_summed.exposure.tolist(), [0.9, 0.8, 0.75, 
                                                             0.75, 0.8, 0.8])

        # different exposure
        exposure = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        bins2 = ExposureBins(self.bins.counts, self.bins.lo_edges, 
                             self.bins.hi_edges, exposure)
        bins_summed = ExposureBins.sum([self.bins, bins2])
        self.assertListEqual(bins_summed.counts.tolist(), [0, 40, 24, 6, 40, 24])
        self.assertListEqual(bins_summed.exposure.tolist(), 
                             [0.95, 0.9, 0.875, 0.875, 0.9, 0.9])
        
        # wrong number of bins
        bins2 = ExposureBins(self.bins.counts[1:], self.bins.lo_edges[1:],
                             self.bins.hi_edges[1:], self.bins.exposure[1:])
        with self.assertRaises(AssertionError):
            bins_summed = ExposureBins.sum([self.bins, bins2])

        # non-matching bin edges
        bins2 = ExposureBins(self.bins.counts, self.bins.lo_edges+0.1,
                             self.bins.hi_edges+0.1, self.bins.exposure)
        with self.assertRaises(AssertionError):
            bins_summed = ExposureBins.sum([self.bins, bins2])

    def test_init_errors(self):
        with self.assertRaises(TypeError):
            ExposureBins(self.bins.counts, self.bins.lo_edges, 
                         self.bins.hi_edges, 1.0)

        with self.assertRaises(ValueError):
            ExposureBins(self.bins.counts, self.bins.lo_edges, 
                         self.bins.hi_edges, self.bins.exposure[1:])


class TestChannelBins(unittest.TestCase):
    
    def setUp(self):
        counts = [20, 50, 17, 3, 0, 3]
        chan_nums = [0, 1, 3, 4, 5, 6]
        exposure = 10.0
        self.bins = ChannelBins.create(counts, chan_nums, exposure)
    
    def test_chan_nums(self):
        self.assertListEqual(self.bins.chan_nums.tolist(), [0, 1, 3, 4, 5, 6])
    
    def test_contiguous_bins(self):
        cont_bins = self.bins.contiguous_bins()
        self.assertEqual(len(cont_bins), 2)
        self.assertTupleEqual(cont_bins[0].range, (0, 1))
        self.assertTupleEqual(cont_bins[1].range, (3, 6))
        
        self.assertEqual(len(cont_bins[0].contiguous_bins()), 1)
    
    def test_rebin(self):
        
        # rebin full range
        rebinned = self.bins.rebin(combine_by_factor, 2)
        self.assertEqual(rebinned.size, 3)
        self.assertListEqual(rebinned.chan_nums.tolist(), [0, 3, 5])
        self.assertListEqual(rebinned.counts.tolist(), [70, 20, 3])
        self.assertListEqual(rebinned.exposure.tolist(), [20.0, 20.0, 20.0])

        # rebin from chan_min through end of range        
        rebinned = self.bins.rebin(combine_by_factor, 2, chan_min=3)
        self.assertEqual(rebinned.size, 4)
        self.assertListEqual(rebinned.chan_nums.tolist(), [0, 1, 3, 5])
        self.assertListEqual(rebinned.counts.tolist(), [20, 50, 20, 3])
        self.assertListEqual(rebinned.exposure.tolist(), [10.0, 10.0, 20.0, 20.0])
        
        # rebin from beginning of range through chan_max
        rebinned = self.bins.rebin(combine_by_factor, 2, chan_max=4)
        self.assertEqual(rebinned.size, 4)
        self.assertListEqual(rebinned.chan_nums.tolist(), [0, 3, 5, 6])
        self.assertListEqual(rebinned.counts.tolist(), [70, 20, 0, 3])
        self.assertListEqual(rebinned.exposure.tolist(), [20.0, 20.0, 10.0, 10.0])

        # rebin middle      
        rebinned = self.bins.rebin(combine_by_factor, 3, chan_min=3, chan_max=5)
        self.assertEqual(rebinned.size, 4)
        self.assertListEqual(rebinned.chan_nums.tolist(), [0, 1, 3, 6])
        self.assertListEqual(rebinned.counts.tolist(), [20, 50, 20, 3])
        self.assertListEqual(rebinned.exposure.tolist(), [10.0, 10.0, 30.0, 10.0])

        # rebin outside range
        rebinned = self.bins.rebin(combine_by_factor, 2, chan_min=7, chan_max=10)
        self.assertEqual(rebinned.size, 6)
        self.assertListEqual(rebinned.chan_nums.tolist(), [0, 1, 3, 4, 5, 6])
        self.assertListEqual(rebinned.counts.tolist(), [20, 50, 17, 3, 0, 3])
        self.assertListEqual(rebinned.exposure.tolist(), [10.0]*6)

    def test_range(self):
        self.assertTupleEqual(self.bins.range, (0, 6))

    def test_slice(self):
        # middle slice       
        bins2 = self.bins.slice(3, 4)
        self.assertTupleEqual(bins2.range, (3, 4))
        
        # slice below lower boundary
        bins2 = self.bins.slice(-1, 3)
        self.assertTupleEqual(bins2.range, (0, 3))

        # slice above upper boundary
        bins2 = self.bins.slice(4, 7)
        self.assertTupleEqual(bins2.range, (4, 6))
        
        # slice covering full range
        bins2 = self.bins.slice(-1, 7)
        self.assertTupleEqual(bins2.range, (0, 6))

        # slice one bin
        bins2 = self.bins.slice(3, 3)
        self.assertTupleEqual(bins2.range, (3, 3))

        # slice fully outside range
        bins2 = self.bins.slice(-2, -1)
        self.assertIsNone(bins2.range)

    def test_merge(self):
        
        counts = [5, 7, 8]
        exposure = [0.9, 0.9, 0.9]
        
        # merge at high end
        chan_nums = [7, 8, 9]
        bins2 = ChannelBins.create(counts, chan_nums, exposure)
        bins_merged = ChannelBins.merge([self.bins, bins2])
        self.assertListEqual(bins_merged.chan_nums.tolist(), [0, 1, 3, 4, 5, 6,
                                                              7, 8, 9])
        
        # flip order
        bins_merged = ChannelBins.merge([bins2, self.bins])
        self.assertListEqual(bins_merged.chan_nums.tolist(), [0, 1, 3, 4, 5, 6,
                                                              7, 8, 9])
        
        # merge at low end
        bins3 = ChannelBins.create([2,3,4], [4,5,6], [0.8, 0.8, 0.8])
        bins_merged = ChannelBins.merge([bins2, bins3])
        self.assertListEqual(bins_merged.chan_nums.tolist(), [4, 5, 6, 7, 8, 9])
        
        # overlapping merge (not allowed)
        bins2 = ChannelBins.create(counts, [5,6,7], exposure)
        with self.assertRaises(ValueError):
            bins_merged = ChannelBins.merge([self.bins, bins2])
        with self.assertRaises(ValueError):
            bins_merged = ChannelBins.merge([bins2, self.bins])

    def test_sum(self):
        
        # same exposure
        bins2 = ChannelBins.create(self.bins.counts, self.bins.chan_nums, 
                                   self.bins.exposure)
        bins_summed = ChannelBins.sum([self.bins, bins2])
        self.assertListEqual(bins_summed.counts.tolist(), [40, 100, 34, 6, 0, 6])
        self.assertListEqual(bins_summed.exposure.tolist(), [10.0]*6)
        
        # different exposure
        exposure = [1.0] * 6
        bins2 = ChannelBins.create(self.bins.counts, self.bins.chan_nums, 
                                   exposure)
        bins_summed = ChannelBins.sum([self.bins, bins2])
        self.assertListEqual(bins_summed.counts.tolist(), [40, 100, 34, 6, 0, 6])
        self.assertListEqual(bins_summed.exposure.tolist(), [5.5]*6)
        
        # wrong number of bins
        bins2 = ChannelBins.create(self.bins.counts[1:], self.bins.chan_nums[1:],
                                   self.bins.exposure[1:])
        with self.assertRaises(AssertionError):
            bins_summed = ChannelBins.sum([self.bins, bins2])

        # non-matching bin edges
        bins2 = ChannelBins.create(self.bins.counts, self.bins.chan_nums + 1,
                                   self.bins.exposure)
        with self.assertRaises(AssertionError):
            bins_summed = ChannelBins.sum([self.bins, bins2])
        
    def test_init_errors(self):
        with self.assertRaises(ValueError):
            ChannelBins.create(self.bins.counts, self.bins.chan_nums[1:], 
                               self.bins.exposure)

        with self.assertRaises(TypeError):
            ChannelBins.create(self.bins.counts, self.bins.chan_nums[0], 
                               self.bins.exposure)
    

class TestTimeBins(unittest.TestCase):
    
    def test_init(self):
        counts = [0, 20, 12, 3]
        exposure = [0.9, 0.8, 0.75, 0.75]
        low = [0.0, 1.0, 2.5, 3.5]
        high = [1.0, 2.0, 3.5, 4.5]
        bins = TimeBins(counts, low, high, exposure)
        self.assertIsInstance(bins, TimeBins)


class TestEnergyBins(unittest.TestCase):
    
    def setUp(self):
        counts = [20, 50, 17, 3, 0]
        low = [10., 50., 180., 320., 500.]    
        high = [50., 180., 320., 500., 1000.]
        exposure = 10.0
        self.bins = EnergyBins(counts, low, high, exposure)
    
    def test_centroids(self):
        vals = [22.36, 94.87, 240., 400., 707.11]
        centroids = self.bins.centroids
        for i in range(4):
            self.assertAlmostEqual(centroids[i], vals[i], places=2)

    def test_rates_per_kev(self):
        vals = [0.050, 0.038, 0.012, 0.002, 0.0]
        r = self.bins.rates_per_kev
        for i in range(4):
            self.assertAlmostEqual(r[i], vals[i], places=3)
    
    def test_rate_uncertainty_per_kev(self):
        vals = [0.0112, 0.0054, 0.0029, 0.0010, 0.0]
        u = self.bins.rate_uncertainty_per_kev
        for i in range(4):
            self.assertAlmostEqual(u[i], vals[i], places=3)

    def test_rebin(self):
        rebinned = self.bins.rebin(combine_by_factor, 2, emax=500.0)
        self.assertEqual(rebinned.size, 3)
        self.assertListEqual(rebinned.hi_edges.tolist(), [180., 500., 1000.])
        self.assertListEqual(rebinned.lo_edges.tolist(), [10., 180., 500.])
        self.assertListEqual(rebinned.counts.tolist(), [70, 20, 0])
        self.assertListEqual(rebinned.exposure.tolist(), [10.0]*3)

    def test_sum(self):
        # same exposure
        bins2 = EnergyBins(self.bins.counts, self.bins.lo_edges, 
                             self.bins.hi_edges, self.bins.exposure)
        bins_summed = EnergyBins.sum([self.bins, bins2])
        self.assertListEqual(bins_summed.counts.tolist(), [40, 100, 34, 6, 0])
        self.assertListEqual(bins_summed.exposure.tolist(), [20.0]*5)
        
        # wrong number of bins
        bins2 = EnergyBins(self.bins.counts[1:], self.bins.lo_edges[1:],
                             self.bins.hi_edges[1:], self.bins.exposure[1:])
        with self.assertRaises(AssertionError):
            bins_summed = EnergyBins.sum([self.bins, bins2])

        # non-matching bin edges
        bins2 = EnergyBins(self.bins.counts, self.bins.lo_edges+0.1,
                             self.bins.hi_edges+0.1, self.bins.exposure)
        with self.assertRaises(AssertionError):
            bins_summed = EnergyBins.sum([self.bins, bins2])


class TestTimeChannelBins(unittest.TestCase):
    
    def setUp(self):
        counts = [[50, 5, 10], [100, 10, 20], [10, 1, 2], [20, 2, 4]]
        tstart = [0.0, 1.0, 3.0, 4.0]
        tstop = [1.0, 2.0, 4.0, 5.0]
        exposure = [1] * 4
        chan_nums = [0, 1, 3]
        self.bins = TimeChannelBins(counts, tstart, tstop, exposure, chan_nums)
    
    def test_chan_nums(self):
        self.assertListEqual(self.bins.chan_nums.tolist(), [0, 1, 3])

    def test_channel_range(self):
        self.assertTupleEqual(self.bins.channel_range, (0, 3))
    
    def test_counts(self):
        self.assertListEqual(self.bins.counts[:,0].tolist(), [50, 100, 10, 20])
        self.assertListEqual(self.bins.counts[:,1].tolist(), [5, 10, 1, 2])
        self.assertListEqual(self.bins.counts[:,2].tolist(), [10, 20, 2, 4])
 
    def test_count_uncertainty(self):
        uncert = self.bins.count_uncertainty
        vals = [7.07, 10., 3.16, 4.47]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,0], vals[i], 2)

        vals = [2.24, 3.16, 1.0, 1.41]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,1], vals[i], 2)

        vals = [3.16, 4.47, 1.41, 2.]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,2], vals[i], 2)
    
    def test_exposure(self):
        self.assertListEqual(self.bins.exposure.tolist(), [1.0]*4)

    def test_num_chans(self):
        self.assertEqual(self.bins.num_chans, 3)

    def test_num_times(self):
        self.assertEqual(self.bins.num_times, 4)

    def test_quality(self):
        self.assertListEqual(self.bins.quality.tolist(), [0]*4)

    def test_rates(self):
        self.assertListEqual(self.bins.rates[:,0].tolist(), [50., 100., 10., 20.])
        self.assertListEqual(self.bins.rates[:,1].tolist(), [5., 10., 1., 2.])
        self.assertListEqual(self.bins.rates[:,2].tolist(), [10., 20., 2., 4.])
        
    def test_rate_uncertainty(self):
        uncert = self.bins.rate_uncertainty
        vals = [7.07, 10., 3.16, 4.47]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,0], vals[i], 2)

        vals = [2.24, 3.16, 1.0, 1.41]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,1], vals[i], 2)

        vals = [3.16, 4.47, 1.41, 2.]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,2], vals[i], 2)

    def test_size(self):
        self.assertTupleEqual(self.bins.size, (4, 3))
    
    def test_time_centroids(self):
        self.assertListEqual(self.bins.time_centroids.tolist(), 
                             [0.5, 1.5, 3.5, 4.5])

    def test_time_range(self):
        self.assertTupleEqual(self.bins.time_range, (0.0, 5.0))
    
    def test_time_widths(self):
        self.assertListEqual(self.bins.time_widths.tolist(), [1.0]*4)
    
    def test_tstart(self):
        self.assertListEqual(self.bins.tstart.tolist(), [0.0, 1.0, 3.0, 4.0])

    def test_tstop(self):
        self.assertListEqual(self.bins.tstop.tolist(), [1.0, 2.0, 4.0, 5.0])

    def test_apply_ebounds(self):
        emin = [10.0, 50.0, 300.0]
        emax = [50.0, 150., 500.0]
        ebounds = Ebounds.from_bounds(emin, emax)
        teb = self.bins.apply_ebounds(ebounds)
        
        assert isinstance(teb, TimeEnergyBins)
        self.assertListEqual(teb.counts.flatten().tolist(), 
                             self.bins.counts.flatten().tolist())
        self.assertListEqual(teb.tstart.tolist(), self.bins.tstart.tolist())
        self.assertListEqual(teb.tstop.tolist(), self.bins.tstop.tolist())
        self.assertListEqual(teb.exposure.tolist(), self.bins.exposure.tolist())
        self.assertListEqual(teb.emin.tolist(), emin)
        self.assertListEqual(teb.emax.tolist(), emax)
        
        with self.assertRaises(TypeError):
            self.bins.apply_ebounds(0)

        with self.assertRaises(ValueError):
            ebounds = Ebounds.from_bounds(emin[1:], emax[1:])
            self.bins.apply_ebounds(ebounds)

    def test_closest_time_edge(self):
        # closest low edge
        self.assertEqual(self.bins.closest_time_edge(0.3, which='low'), 0.0)
        # closest high edge
        self.assertEqual(self.bins.closest_time_edge(0.3, which='high'), 1.0)
        # closest edge
        self.assertEqual(self.bins.closest_time_edge(0.3, which='either'), 0.0)

    def test_contiguous_channel_bins(self):
        cont_bins = self.bins.contiguous_channel_bins()
        self.assertEqual(len(cont_bins), 2)
        self.assertTupleEqual(cont_bins[0].channel_range, (0, 1))
        self.assertTupleEqual(cont_bins[1].channel_range, (3, 3))
        
        self.assertEqual(len(cont_bins[0].contiguous_channel_bins()), 1)

    def test_contiguous_time_bins(self):
        cont_bins = self.bins.contiguous_time_bins()
        self.assertEqual(len(cont_bins), 2)
        self.assertTupleEqual(cont_bins[0].time_range, (0., 2.))
        self.assertTupleEqual(cont_bins[1].time_range, (3., 5.))
        
        self.assertEqual(len(cont_bins[0].contiguous_time_bins()), 1)
    
    def test_get_exposure(self):
        # full range
        self.assertEqual(self.bins.get_exposure(), 4.0)
        
        # one time range
        self.assertEqual(self.bins.get_exposure(time_ranges=(2.0, 4.0)), 2.0)
        
        # two time ranges
        self.assertEqual(self.bins.get_exposure(time_ranges=[(0.0, 1.0), 
                                                             (3.0, 4.0)]), 2.0)

        # offset edges, no scaling
        self.assertEqual(self.bins.get_exposure(time_ranges=(2.5, 4.0)), 2.0)

        # offset edges, with scaling
        self.assertEqual(self.bins.get_exposure(time_ranges=(2.5, 4.0), 
                                                scale=True), 1.5)

    def test_integrate_channels(self):
        # full range
        bins = self.bins.integrate_channels()
        self.assertIsInstance(bins, TimeBins)
        self.assertListEqual(bins.counts.tolist(), [65, 130, 13, 26])
        self.assertListEqual(bins.exposure.tolist(), [1.0]*4)
        self.assertListEqual(bins.lo_edges.tolist(), [0.0, 1.0, 3.0, 4.0])
        self.assertListEqual(bins.hi_edges.tolist(), [1.0, 2.0, 4.0, 5.0])

        # set emin
        bins = self.bins.integrate_channels(chan_min=1)
        self.assertListEqual(bins.counts.tolist(), [15, 30, 3, 6])
        self.assertListEqual(bins.exposure.tolist(), [1.0]*4)
        self.assertListEqual(bins.lo_edges.tolist(), [0.0, 1.0, 3.0, 4.0])
        self.assertListEqual(bins.hi_edges.tolist(), [1.0, 2.0, 4.0, 5.0])

        # set emax
        bins = self.bins.integrate_channels(chan_max=1)
        self.assertListEqual(bins.counts.tolist(), [55, 110, 11, 22])
        self.assertListEqual(bins.exposure.tolist(), [1.0]*4)
        self.assertListEqual(bins.lo_edges.tolist(), [0.0, 1.0, 3.0, 4.0])
        self.assertListEqual(bins.hi_edges.tolist(), [1.0, 2.0, 4.0, 5.0])

        # set both emin and emax
        bins = self.bins.integrate_channels(chan_min=1, chan_max=1)
        self.assertListEqual(bins.counts.tolist(), [5, 10, 1, 2])
        self.assertListEqual(bins.exposure.tolist(), [1.0]*4)
        self.assertListEqual(bins.lo_edges.tolist(), [0.0, 1.0, 3.0, 4.0])
        self.assertListEqual(bins.hi_edges.tolist(), [1.0, 2.0, 4.0, 5.0])
    
    def test_integrate_time(self):
        # full range
        bins = self.bins.integrate_time()
        self.assertIsInstance(bins, ChannelBins)
        self.assertListEqual(bins.counts.tolist(), [180, 18, 36])
        self.assertListEqual(bins.exposure.tolist(), [4.0]*3)
        self.assertListEqual(bins.chan_nums.tolist(), [0, 1, 3])
        
        # set tstart
        bins = self.bins.integrate_time(tstart=1.5)
        self.assertListEqual(bins.counts.tolist(), [130, 13, 26])
        self.assertListEqual(bins.exposure.tolist(), [3.0]*3)
        self.assertListEqual(bins.chan_nums.tolist(), [0, 1, 3])

        # set tstop
        bins = self.bins.integrate_time(tstop=3.5)
        self.assertListEqual(bins.counts.tolist(), [160, 16, 32])
        self.assertListEqual(bins.exposure.tolist(), [3.0]*3)
        self.assertListEqual(bins.chan_nums.tolist(), [0, 1, 3])

        # set both tstart amd tstop
        bins = self.bins.integrate_time(tstart=1.5, tstop=3.5)
        self.assertListEqual(bins.counts.tolist(), [110, 11, 22])
        self.assertListEqual(bins.exposure.tolist(), [2.0]*3)
        self.assertListEqual(bins.chan_nums.tolist(), [0, 1, 3])

    def test_rebin_channels(self):
        
        # rebin full range
        rebinned = self.bins.rebin_channels(combine_by_factor, 2)
        self.assertEqual(rebinned.num_chans, 1)
        self.assertListEqual(rebinned.chan_nums.tolist(), [0])
        counts = rebinned.counts.tolist()
        vals = [55, 110, 11, 22]
        for i in range(4):
            self.assertEqual(counts[i][0], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
        
        # rebin from emin through end of range        
        rebinned = self.bins.rebin_channels(combine_by_factor, 2, chan_min=0)
        self.assertEqual(rebinned.num_chans, 1)
        self.assertListEqual(rebinned.chan_nums.tolist(), [0])
        counts = rebinned.counts.tolist()
        vals = [55, 110, 11, 22]
        for i in range(4):
            self.assertEqual(counts[i][0], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
        
        # rebin from beginning of range through emax       
        rebinned = self.bins.rebin_channels(combine_by_factor, 2, chan_max=3)
        self.assertEqual(rebinned.num_chans, 1)
        self.assertListEqual(rebinned.chan_nums.tolist(), [0])
        counts = rebinned.counts.tolist()
        vals = [55, 110, 11, 22]
        for i in range(4):
            self.assertEqual(counts[i][0], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
        
        # rebin middle      
        rebinned = self.bins.rebin_channels(combine_by_factor, 2, chan_min=0, 
                                            chan_max=1)
        self.assertEqual(rebinned.num_chans, 2)
        self.assertListEqual(rebinned.chan_nums.tolist(), [0, 3])
        counts = rebinned.counts.tolist()
        vals = [[55, 10], [110, 20], [11, 2], [22, 4]]
        for i in range(4):
            self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
 
        # rebin outside range
        rebinned = self.bins.rebin_channels(combine_by_factor, 2, chan_min=4, 
                                            chan_max=6)
        self.assertEqual(rebinned.num_chans, 3)
        self.assertListEqual(rebinned.chan_nums.tolist(), [0, 1, 3])
        counts = rebinned.counts.tolist()
        vals = [[50, 5, 10], [100, 10, 20], [10, 1, 2], [20, 2, 4]]
        for i in range(4):
            self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
   
    def test_rebin_time(self):
        # rebin full range
        rebinned = self.bins.rebin_time(combine_by_factor, 2)
        self.assertEqual(rebinned.num_times, 2)
        self.assertListEqual(rebinned.tstart.tolist(), [0.0, 3.0])
        self.assertListEqual(rebinned.tstop.tolist(), [2.0, 5.0])
        counts = rebinned.counts.tolist()
        vals = [[150, 15, 30], [30, 3, 6]]
        for i in range(2):
            self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [2.0, 2.0])
        
        # rebin from tstart through end of range        
        rebinned = self.bins.rebin_time(combine_by_factor, 2, tstart=0.5)
        self.assertEqual(rebinned.num_times, 2)
        self.assertListEqual(rebinned.tstart.tolist(), [0.0, 3.0])
        self.assertListEqual(rebinned.tstop.tolist(), [2.0, 5.0])
        counts = rebinned.counts.tolist()
        vals = [[150, 15, 30], [30, 3, 6]]
        for i in range(2):
             self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [2.0, 2.0])
        
        # rebin from beginning of range through tstop       
        rebinned = self.bins.rebin_time(combine_by_factor, 2, tstop=4.5)
        self.assertEqual(rebinned.num_times, 2)
        self.assertListEqual(rebinned.tstart.tolist(), [0.0, 3.0])
        self.assertListEqual(rebinned.tstop.tolist(), [2.0, 5.0])
        counts = rebinned.counts.tolist()
        vals = [[150, 15, 30], [30, 3, 6]]
        for i in range(2):
             self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [2.0, 2.0])
        
        # rebin middle      
        rebinned = self.bins.rebin_time(combine_by_factor, 2, tstart=0.5, 
                                          tstop=1.5)
        self.assertEqual(rebinned.num_times, 3)
        self.assertListEqual(rebinned.tstart.tolist(), [0.0, 3.0, 4.0])
        self.assertListEqual(rebinned.tstop.tolist(), [2.0, 4.0, 5.0])
        counts = rebinned.counts.tolist()
        vals = [[150, 15, 30], [10, 1, 2], [20, 2, 4]]
        for i in range(3):
            self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [2.0, 1.0, 1.0])

 
        # rebin outside range
        rebinned = self.bins.rebin_time(combine_by_factor, 2, tstart=10.0, 
                                        tstop=20.0)
        self.assertEqual(rebinned.num_times, 4)
        self.assertListEqual(rebinned.tstart.tolist(), [0.0, 1.0, 3.0, 4.0])
        self.assertListEqual(rebinned.tstop.tolist(), [1.0, 2.0, 4.0, 5.0])
        counts = rebinned.counts.tolist()
        vals = [[50, 5, 10], [100, 10, 20], [10, 1, 2], [20, 2, 4]]
        for i in range(4):
            self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
    
    def test_slice_channels(self):
        # middle slice       
        bins2 = self.bins.slice_channels(0, 1)
        self.assertTupleEqual(bins2.channel_range, (0, 1))
        
        # slice below lower boundary
        bins2 = self.bins.slice_channels(-1, 1)
        self.assertTupleEqual(bins2.channel_range, (0, 1))

        # slice above upper boundary
        bins2 = self.bins.slice_channels(1, 5)
        self.assertTupleEqual(bins2.channel_range, (1, 3))
        
        # slice covering full range
        bins2 = self.bins.slice_channels(-1, 5)
        self.assertTupleEqual(bins2.channel_range, (0, 3))

        # slice one bin
        bins2 = self.bins.slice_channels(1, 1)
        self.assertTupleEqual(bins2.channel_range, (1, 1))

        # slice fully outside range
        bins2 = self.bins.slice_channels(5, 10)
        self.assertIsNone(bins2.channel_range)

    def test_slice_time(self):
        # middle slice       
        bins2 = self.bins.slice_time(1.5, 3.5)
        self.assertTupleEqual(bins2.time_range, (1.0, 4.0))
        
        # slice below lower boundary
        bins2 = self.bins.slice_time(-1.0, 3.5)
        self.assertTupleEqual(bins2.time_range, (0.0, 4.0))

        # slice above upper boundary
        bins2 = self.bins.slice_time(1.5, 10.0)
        self.assertTupleEqual(bins2.time_range, (1.0, 5.0))
        
        # slice covering full range
        bins2 = self.bins.slice_time(-1.0, 10.0)
        self.assertTupleEqual(bins2.time_range, (0.0, 5.0))

        # slice one bin
        bins2 = self.bins.slice_time(1.5, 1.5)
        self.assertTupleEqual(bins2.time_range, (1.0, 2.0))

        # slice fully outside range
        bins2 = self.bins.slice_time(10.0, 20.0)
        self.assertIsNone(bins2.time_range)

    def test_merge_channels(self):
        
        counts = [[10, 20], [10, 20], [10, 20], [10, 20]]        
        # merge at high end
        chan_nums = [4, 5]
        bins2 = TimeChannelBins(counts, self.bins.tstart, self.bins.tstop, 
                                self.bins.exposure, chan_nums)
        bins_merged = TimeChannelBins.merge_channels([self.bins, bins2])
        self.assertListEqual(bins_merged.chan_nums.tolist(), [0, 1, 3, 4, 5])
        
        # flip order
        bins_merged = TimeChannelBins.merge_channels([bins2, self.bins])
        self.assertListEqual(bins_merged.chan_nums.tolist(), [0, 1, 3, 4, 5])
        
        # merge at low end
        bins3 = TimeChannelBins(counts, self.bins.tstart, self.bins.tstop, 
                                self.bins.exposure, [2, 3])
        bins_merged = TimeChannelBins.merge_channels([bins2, bins3])
        self.assertListEqual(bins_merged.chan_nums.tolist(), [2, 3, 4, 5])
                
        # overlapping merge (not allowed)
        chan_nums = [2, 3]
        bins2 = TimeChannelBins(counts, self.bins.tstart, self.bins.tstop, 
                                self.bins.exposure, chan_nums)
        with self.assertRaises(ValueError):
            bins_merged = TimeChannelBins.merge_channels([self.bins, bins2])
        with self.assertRaises(ValueError):
            bins_merged = TimeChannelBins.merge_channels([bins2, self.bins])

    def test_merge_time(self):
        
        counts = [[10, 20, 30], [10, 20, 30]] 
        exposure = [1.0, 1.0]       
        # merge at high end
        tstart = [5.0, 6.0]
        tstop = [6.0, 7.0]
        bins2 = TimeChannelBins(counts, tstart, tstop, exposure, 
                                self.bins.chan_nums)
        bins_merged = TimeChannelBins.merge_time([self.bins, bins2])
        self.assertListEqual(bins_merged.tstop.tolist(), [1.0, 2.0, 4.0, 5.0, 
                                                          6.0, 7.0])
        self.assertListEqual(bins_merged.tstart.tolist(), [0.0, 1.0, 3.0, 4.0,
                                                           5.0, 6.0])
        # flip order
        bins_merged = TimeChannelBins.merge_time([bins2, self.bins])
        self.assertListEqual(bins_merged.tstop.tolist(), [1.0, 2.0, 4.0, 5.0, 
                                                          6.0, 7.0])
        self.assertListEqual(bins_merged.tstart.tolist(), [0.0, 1.0, 3.0, 4.0,
                                                           5.0, 6.0])
        # merge at low end
        tstart = [-2.0, -1.0]
        tstop = [-1.0, 0.0]
        bins2 = TimeChannelBins(counts, tstart, tstop, exposure, 
                                self.bins.chan_nums)
        bins_merged = TimeChannelBins.merge_time([self.bins, bins2])
        self.assertListEqual(bins_merged.tstop.tolist(), [-1.0, 0.0, 1.0, 2.0, 
                                                          4.0, 5.0])
        self.assertListEqual(bins_merged.tstart.tolist(), [-2.0, -1.0, 0.0, 1.0,
                                                           3.0, 4.0])
                
        # overlapping merge (not allowed)
        tstart = [4.0, 5.0]
        tstop = [5.0, 6.0]
        bins2 = TimeChannelBins(counts, tstart, tstop, exposure, 
                                self.bins.chan_nums)
        with self.assertRaises(ValueError):
            bins_merged = TimeChannelBins.merge_time([self.bins, bins2])
        with self.assertRaises(ValueError):
            bins_merged = TimeChannelBins.merge_time([bins2, self.bins])

    def test_init_errors(self):
        with self.assertRaises(TypeError):
            TimeChannelBins(0, self.bins.tstart, self.bins.tstop, 
                            self.bins.exposure, self.bins.chan_nums)

        with self.assertRaises(TypeError):
            TimeEnergyBins(self.bins.counts[:,0], self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure, 
                           self.bins.chan_nums)
        
        with self.assertRaises(TypeError):
            TimeChannelBins(self.bins.counts, self.bins.tstart[0], 
                            self.bins.tstop, self.bins.exposure, 
                            self.bins.chan_nums)

        with self.assertRaises(TypeError):
            TimeChannelBins(self.bins.counts, self.bins.tstart, 
                            self.bins.tstop[0], self.bins.exposure, 
                            self.bins.chan_nums)

        with self.assertRaises(TypeError):
            TimeChannelBins(self.bins.counts, self.bins.tstart, 
                            self.bins.tstop, self.bins.exposure[0], 
                            self.bins.chan_nums)

        with self.assertRaises(TypeError):
            TimeChannelBins(self.bins.counts, self.bins.tstart, 
                            self.bins.tstop, self.bins.exposure, 
                            self.bins.chan_nums[0])

        with self.assertRaises(ValueError):
            TimeChannelBins(self.bins.counts, self.bins.tstart[1:], 
                            self.bins.tstop, self.bins.exposure, 
                            self.bins.chan_nums)

        with self.assertRaises(ValueError):
            TimeChannelBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure[1:], 
                           self.bins.chan_nums)

        with self.assertRaises(ValueError):
            TimeChannelBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure[1:], 
                           self.bins.chan_nums[1:])

        with self.assertRaises(ValueError):
            TimeChannelBins(self.bins.counts.T, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure[1:], 
                           self.bins.chan_nums)

        with self.assertRaises(TypeError):
            TimeChannelBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure, 
                           self.bins.chan_nums, quality=0)

        with self.assertRaises(ValueError):
            TimeChannelBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure, 
                           self.bins.chan_nums, quality=[0, 0])
 

class TestTimeEnergyBins(unittest.TestCase):
    
    def setUp(self):
        counts = [[50, 5, 10], [100, 10, 20], [10, 1, 2], [20, 2, 4]]
        tstart = [0.0, 1.0, 3.0, 4.0]
        tstop = [1.0, 2.0, 4.0, 5.0]
        exposure = [1] * 4
        emin = [10.0, 50.0, 300.0]
        emax = [50.0, 150., 500.0]
        self.bins = TimeEnergyBins(counts, tstart, tstop, exposure, emin, emax)
    
    def test_chan_widths(self):
        self.assertListEqual(self.bins.chan_widths.tolist(), [40., 100., 200.])
    
    def test_counts(self):
        self.assertListEqual(self.bins.counts[:,0].tolist(), [50, 100, 10, 20])
        self.assertListEqual(self.bins.counts[:,1].tolist(), [5, 10, 1, 2])
        self.assertListEqual(self.bins.counts[:,2].tolist(), [10, 20, 2, 4])
 
    def test_count_uncertainty(self):
        uncert = self.bins.count_uncertainty
        vals = [7.07, 10., 3.16, 4.47]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,0], vals[i], 2)

        vals = [2.24, 3.16, 1.0, 1.41]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,1], vals[i], 2)

        vals = [3.16, 4.47, 1.41, 2.]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,2], vals[i], 2)

    def test_emax(self):
        self.assertListEqual(self.bins.emax.tolist(), [50., 150., 500.])
    
    def test_emin(self):
        self.assertListEqual(self.bins.emin.tolist(), [10., 50., 300.])

    def test_energy_centroids(self):
        vals = [22.36, 86.60, 387.30]
        centroids = self.bins.energy_centroids
        for i in range(3):
            self.assertAlmostEqual(centroids[i], vals[i], places=2)
    
    def test_energy_range(self):
        self.assertTupleEqual(self.bins.energy_range, (10.0, 500.0))

    def test_exposure(self):
        self.assertListEqual(self.bins.exposure.tolist(), [1.0]*4)

    def test_num_chans(self):
        self.assertEqual(self.bins.num_chans, 3)

    def test_num_times(self):
        self.assertEqual(self.bins.num_times, 4)

    def test_quality(self):
        self.assertListEqual(self.bins.quality.tolist(), [0]*4)

    def test_rates(self):
        self.assertListEqual(self.bins.rates[:,0].tolist(), [50., 100., 10., 20.])
        self.assertListEqual(self.bins.rates[:,1].tolist(), [5., 10., 1., 2.])
        self.assertListEqual(self.bins.rates[:,2].tolist(), [10., 20., 2., 4.])

    def test_rates_per_kev(self):
        r = self.bins.rates_per_kev
        vals = [1.25, 2.5, 0.25, 0.5]
        for i in range(4):
            self.assertAlmostEqual(r[i,0], vals[i], places=2)
        vals = [0.05, 0.1, 0.01, 0.02]
        for i in range(4):
            self.assertAlmostEqual(r[i,1], vals[i], places=2)
        vals = [0.05, 0.1, 0.01, 0.02]
        for i in range(4):
            self.assertAlmostEqual(r[i,2], vals[i], places=2)
        
    def test_rate_uncertainty(self):
        uncert = self.bins.rate_uncertainty
        vals = [7.07, 10., 3.16, 4.47]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,0], vals[i], 2)

        vals = [2.24, 3.16, 1.0, 1.41]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,1], vals[i], 2)

        vals = [3.16, 4.47, 1.41, 2.]
        for i in range(4):
            self.assertAlmostEqual(uncert[i,2], vals[i], 2)

    def test_rate_uncertainty_per_kev(self):
        uncert = self.bins.rate_uncertainty_per_kev
        vals = [0.177, 0.25, 0.079, 0.112] 
        for i in range(4):
            self.assertAlmostEqual(uncert[i,0], vals[i], 3)

        vals = [0.022, 0.032, 0.01, 0.014] 
        for i in range(4):
            self.assertAlmostEqual(uncert[i,1], vals[i], 3)

        vals = [0.016, 0.022, 0.007, 0.01] 
        for i in range(4):
            self.assertAlmostEqual(uncert[i,2], vals[i], 3)

    def test_size(self):
        self.assertTupleEqual(self.bins.size, (4, 3))
    
    def test_time_centroids(self):
        self.assertListEqual(self.bins.time_centroids.tolist(), 
                             [0.5, 1.5, 3.5, 4.5])

    def test_time_range(self):
        self.assertTupleEqual(self.bins.time_range, (0.0, 5.0))
    
    def test_time_widths(self):
        self.assertListEqual(self.bins.time_widths.tolist(), [1.0]*4)
    
    def test_tstart(self):
        self.assertListEqual(self.bins.tstart.tolist(), [0.0, 1.0, 3.0, 4.0])

    def test_tstop(self):
        self.assertListEqual(self.bins.tstop.tolist(), [1.0, 2.0, 4.0, 5.0])

    def test_closest_energy_edge(self):
        # closest low edge
        self.assertEqual(self.bins.closest_energy_edge(20.0, which='low'), 10.0)
        # closest high edge
        self.assertEqual(self.bins.closest_energy_edge(20.0, which='high'), 50.0)
        # closest edge
        self.assertEqual(self.bins.closest_energy_edge(20.0, which='either'), 10.0)

    def test_closest_time_edge(self):
        # closest low edge
        self.assertEqual(self.bins.closest_time_edge(0.3, which='low'), 0.0)
        # closest high edge
        self.assertEqual(self.bins.closest_time_edge(0.3, which='high'), 1.0)
        # closest edge
        self.assertEqual(self.bins.closest_time_edge(0.3, which='either'), 0.0)

    def test_contiguous_energy_bins(self):
        cont_bins = self.bins.contiguous_energy_bins()
        self.assertEqual(len(cont_bins), 2)
        self.assertTupleEqual(cont_bins[0].energy_range, (10., 150.))
        self.assertTupleEqual(cont_bins[1].energy_range, (300., 500.))
        
        self.assertEqual(len(cont_bins[0].contiguous_energy_bins()), 1)

    def test_contiguous_time_bins(self):
        cont_bins = self.bins.contiguous_time_bins()
        self.assertEqual(len(cont_bins), 2)
        self.assertTupleEqual(cont_bins[0].time_range, (0., 2.))
        self.assertTupleEqual(cont_bins[1].time_range, (3., 5.))
        
        self.assertEqual(len(cont_bins[0].contiguous_time_bins()), 1)
    
    def test_get_exposure(self):
        # full range
        self.assertEqual(self.bins.get_exposure(), 4.0)
        
        # one time range
        self.assertEqual(self.bins.get_exposure(time_ranges=(2.0, 4.0)), 2.0)
        
        # two time ranges
        self.assertEqual(self.bins.get_exposure(time_ranges=[(0.0, 1.0), 
                                                             (3.0, 4.0)]), 2.0)

        # offset edges, no scaling
        self.assertEqual(self.bins.get_exposure(time_ranges=(2.5, 4.0)), 2.0)

        # offset edges, with scaling
        self.assertEqual(self.bins.get_exposure(time_ranges=(2.5, 4.0), 
                                                scale=True), 1.5)

    def test_integrate_energy(self):
        # full range
        bins = self.bins.integrate_energy()
        self.assertIsInstance(bins, TimeBins)
        self.assertListEqual(bins.counts.tolist(), [65, 130, 13, 26])
        self.assertListEqual(bins.exposure.tolist(), [1.0]*4)
        self.assertListEqual(bins.lo_edges.tolist(), [0.0, 1.0, 3.0, 4.0])
        self.assertListEqual(bins.hi_edges.tolist(), [1.0, 2.0, 4.0, 5.0])

        # set emin
        bins = self.bins.integrate_energy(emin=100.0)
        self.assertListEqual(bins.counts.tolist(), [15, 30, 3, 6])
        self.assertListEqual(bins.exposure.tolist(), [1.0]*4)
        self.assertListEqual(bins.lo_edges.tolist(), [0.0, 1.0, 3.0, 4.0])
        self.assertListEqual(bins.hi_edges.tolist(), [1.0, 2.0, 4.0, 5.0])

        # set emax
        bins = self.bins.integrate_energy(emax=100.0)
        self.assertListEqual(bins.counts.tolist(), [55, 110, 11, 22])
        self.assertListEqual(bins.exposure.tolist(), [1.0]*4)
        self.assertListEqual(bins.lo_edges.tolist(), [0.0, 1.0, 3.0, 4.0])
        self.assertListEqual(bins.hi_edges.tolist(), [1.0, 2.0, 4.0, 5.0])

        # set both emin and emax
        bins = self.bins.integrate_energy(emin=75., emax=125.)
        self.assertListEqual(bins.counts.tolist(), [5, 10, 1, 2])
        self.assertListEqual(bins.exposure.tolist(), [1.0]*4)
        self.assertListEqual(bins.lo_edges.tolist(), [0.0, 1.0, 3.0, 4.0])
        self.assertListEqual(bins.hi_edges.tolist(), [1.0, 2.0, 4.0, 5.0])

    def test_integrate_time(self):
        # full range
        bins = self.bins.integrate_time()
        self.assertIsInstance(bins, EnergyBins)
        self.assertListEqual(bins.counts.tolist(), [180, 18, 36])
        self.assertListEqual(bins.exposure.tolist(), [4.0]*3)
        self.assertListEqual(bins.lo_edges.tolist(), [10., 50., 300.])
        self.assertListEqual(bins.hi_edges.tolist(), [50., 150., 500.])

        # set tstart
        bins = self.bins.integrate_time(tstart=1.5)
        self.assertListEqual(bins.counts.tolist(), [130, 13, 26])
        self.assertListEqual(bins.exposure.tolist(), [3.0]*3)
        self.assertListEqual(bins.lo_edges.tolist(), [10., 50., 300.])
        self.assertListEqual(bins.hi_edges.tolist(), [50., 150., 500.])

        # set tstop
        bins = self.bins.integrate_time(tstop=3.5)
        self.assertListEqual(bins.counts.tolist(), [160, 16, 32])
        self.assertListEqual(bins.exposure.tolist(), [3.0]*3)
        self.assertListEqual(bins.lo_edges.tolist(), [10., 50., 300.])
        self.assertListEqual(bins.hi_edges.tolist(), [50., 150., 500.])

        # set both tstart amd tstop
        bins = self.bins.integrate_time(tstart=1.5, tstop=3.5)
        self.assertListEqual(bins.counts.tolist(), [110, 11, 22])
        self.assertListEqual(bins.exposure.tolist(), [2.0]*3)
        self.assertListEqual(bins.lo_edges.tolist(), [10., 50., 300.])
        self.assertListEqual(bins.hi_edges.tolist(), [50., 150., 500.])

    def test_rebin_energy(self):
        
        # rebin full range
        rebinned = self.bins.rebin_energy(combine_by_factor, 2)
        self.assertEqual(rebinned.num_chans, 1)
        self.assertListEqual(rebinned.emax.tolist(), [150.])
        self.assertListEqual(rebinned.emin.tolist(), [10.])
        counts = rebinned.counts.tolist()
        vals = [55, 110, 11, 22]
        for i in range(4):
            self.assertEqual(counts[i][0], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
        
        # rebin from emin through end of range        
        rebinned = self.bins.rebin_energy(combine_by_factor, 2, emin=20.)
        self.assertEqual(rebinned.num_chans, 1)
        self.assertListEqual(rebinned.emax.tolist(), [150.])
        self.assertListEqual(rebinned.emin.tolist(), [10.])
        counts = rebinned.counts.tolist()
        vals = [55, 110, 11, 22]
        for i in range(4):
            self.assertEqual(counts[i][0], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
        
        # rebin from beginning of range through emax       
        rebinned = self.bins.rebin_energy(combine_by_factor, 2, emax=350.0)
        self.assertEqual(rebinned.num_chans, 1)
        self.assertListEqual(rebinned.emax.tolist(), [150.])
        self.assertListEqual(rebinned.emin.tolist(), [10.])
        counts = rebinned.counts.tolist()
        vals = [55, 110, 11, 22]
        for i in range(4):
            self.assertEqual(counts[i][0], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
        
        # rebin middle      
        rebinned = self.bins.rebin_energy(combine_by_factor, 2, emin=20.0, 
                                          emax=125.0)
        self.assertEqual(rebinned.num_chans, 2)
        self.assertListEqual(rebinned.emax.tolist(), [150., 500.])
        self.assertListEqual(rebinned.emin.tolist(), [10., 300.])
        counts = rebinned.counts.tolist()
        vals = [[55, 10], [110, 20], [11, 2], [22, 4]]
        for i in range(4):
            self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
 
        # rebin outside range
        rebinned = self.bins.rebin_energy(combine_by_factor, 2, emin=1000.0, 
                                   emax=2000.)
        self.assertEqual(rebinned.num_chans, 3)
        self.assertListEqual(rebinned.emax.tolist(), [50.0, 150., 500.0])
        self.assertListEqual(rebinned.emin.tolist(), [10.0, 50.0, 300.0])
        counts = rebinned.counts.tolist()
        vals = [[50, 5, 10], [100, 10, 20], [10, 1, 2], [20, 2, 4]]
        for i in range(4):
            self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
   
    def test_rebin_time(self):
        # rebin full range
        rebinned = self.bins.rebin_time(combine_by_factor, 2)
        self.assertEqual(rebinned.num_times, 2)
        self.assertListEqual(rebinned.tstart.tolist(), [0.0, 3.0])
        self.assertListEqual(rebinned.tstop.tolist(), [2.0, 5.0])
        counts = rebinned.counts.tolist()
        vals = [[150, 15, 30], [30, 3, 6]]
        for i in range(2):
            self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [2.0, 2.0])
        
        # rebin from emin through end of range        
        rebinned = self.bins.rebin_time(combine_by_factor, 2, tstart=0.5)
        self.assertEqual(rebinned.num_times, 2)
        self.assertListEqual(rebinned.tstart.tolist(), [0.0, 3.0])
        self.assertListEqual(rebinned.tstop.tolist(), [2.0, 5.0])
        counts = rebinned.counts.tolist()
        vals = [[150, 15, 30], [30, 3, 6]]
        for i in range(2):
             self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [2.0, 2.0])
        
        # rebin from beginning of range through emax       
        rebinned = self.bins.rebin_time(combine_by_factor, 2, tstop=4.5)
        self.assertEqual(rebinned.num_times, 2)
        self.assertListEqual(rebinned.tstart.tolist(), [0.0, 3.0])
        self.assertListEqual(rebinned.tstop.tolist(), [2.0, 5.0])
        counts = rebinned.counts.tolist()
        vals = [[150, 15, 30], [30, 3, 6]]
        for i in range(2):
             self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [2.0, 2.0])
        
        # rebin middle      
        rebinned = self.bins.rebin_time(combine_by_factor, 2, tstart=0.5, 
                                          tstop=1.5)
        self.assertEqual(rebinned.num_times, 3)
        self.assertListEqual(rebinned.tstart.tolist(), [0.0, 3.0, 4.0])
        self.assertListEqual(rebinned.tstop.tolist(), [2.0, 4.0, 5.0])
        counts = rebinned.counts.tolist()
        vals = [[150, 15, 30], [10, 1, 2], [20, 2, 4]]
        for i in range(3):
            self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [2.0, 1.0, 1.0])

 
        # rebin outside range
        rebinned = self.bins.rebin_time(combine_by_factor, 2, tstart=10.0, 
                                        tstop=20.0)
        self.assertEqual(rebinned.num_times, 4)
        self.assertListEqual(rebinned.tstart.tolist(), [0.0, 1.0, 3.0, 4.0])
        self.assertListEqual(rebinned.tstop.tolist(), [1.0, 2.0, 4.0, 5.0])
        counts = rebinned.counts.tolist()
        vals = [[50, 5, 10], [100, 10, 20], [10, 1, 2], [20, 2, 4]]
        for i in range(4):
            self.assertListEqual(counts[i], vals[i])
        self.assertListEqual(rebinned.exposure.tolist(), [1.0]*4)
    
    def test_slice_energy(self):
        # middle slice       
        bins2 = self.bins.slice_energy(20.0, 125.0)
        self.assertTupleEqual(bins2.energy_range, (10.0, 150.0))
        
        # slice below lower boundary
        bins2 = self.bins.slice_energy(1.0, 125.0)
        self.assertTupleEqual(bins2.energy_range, (10.0, 150.0))

        # slice above upper boundary
        bins2 = self.bins.slice_energy(125.0, 1000.0)
        self.assertTupleEqual(bins2.energy_range, (50.0, 500.0))
        
        # slice covering full range
        bins2 = self.bins.slice_energy(1.0, 1000.0)
        self.assertTupleEqual(bins2.energy_range, (10.0, 500.0))

        # slice one bin
        bins2 = self.bins.slice_energy(70.0, 70.0)
        self.assertTupleEqual(bins2.energy_range, (50.0, 150.0))

        # slice fully outside range
        bins2 = self.bins.slice_energy(1000.0, 2000.0)
        self.assertIsNone(bins2.energy_range)

    def test_slice_time(self):
        # middle slice       
        bins2 = self.bins.slice_time(1.5, 3.5)
        self.assertTupleEqual(bins2.time_range, (1.0, 4.0))
        
        # slice below lower boundary
        bins2 = self.bins.slice_time(-1.0, 3.5)
        self.assertTupleEqual(bins2.time_range, (0.0, 4.0))

        # slice above upper boundary
        bins2 = self.bins.slice_time(1.5, 10.0)
        self.assertTupleEqual(bins2.time_range, (1.0, 5.0))
        
        # slice covering full range
        bins2 = self.bins.slice_time(-1.0, 10.0)
        self.assertTupleEqual(bins2.time_range, (0.0, 5.0))

        # slice one bin
        bins2 = self.bins.slice_time(1.5, 1.5)
        self.assertTupleEqual(bins2.time_range, (1.0, 2.0))

        # slice fully outside range
        bins2 = self.bins.slice_time(10.0, 20.0)
        self.assertIsNone(bins2.time_range)

    def test_merge_energy(self):
        
        counts = [[10, 20], [10, 20], [10, 20], [10, 20]]        
        # merge at high end
        emin = [500.0, 1000.0]
        emax = [1000.0, 2000.0]
        bins2 = TimeEnergyBins(counts, self.bins.tstart, self.bins.tstop, 
                               self.bins.exposure, emin, emax)
        bins_merged = TimeEnergyBins.merge_energy([self.bins, bins2])
        self.assertListEqual(bins_merged.emax.tolist(), [50., 150., 500., 1000.,
                                                         2000.])
        self.assertListEqual(bins_merged.emin.tolist(), [10., 50., 300., 500.,
                                                         1000.])
        
        # flip order
        bins_merged = TimeEnergyBins.merge_energy([bins2, self.bins])
        self.assertListEqual(bins_merged.emax.tolist(), [50., 150., 500., 1000.,
                                                         2000.])
        self.assertListEqual(bins_merged.emin.tolist(), [10., 50., 300., 500.,
                                                         1000.])
        # merge at low end
        emin = [1.0, 3.0]
        emax = [3.0, 10.0]
        bins2 = TimeEnergyBins(counts, self.bins.tstart, self.bins.tstop, 
                               self.bins.exposure, emin, emax)
        bins_merged = TimeEnergyBins.merge_energy([self.bins, bins2])
        self.assertListEqual(bins_merged.emax.tolist(),[3., 10., 50., 150., 500.,])
        self.assertListEqual(bins_merged.emin.tolist(), [1., 3., 10., 50., 300.])
                
        # overlapping merge (not allowed)
        emin = [200.0, 700.0]
        emax = [700.0, 1500.0]
        bins2 = TimeEnergyBins(counts, self.bins.tstart, self.bins.tstop, 
                               self.bins.exposure, emin, emax)
        with self.assertRaises(ValueError):
            bins_merged = TimeEnergyBins.merge_energy([self.bins, bins2])
        with self.assertRaises(ValueError):
            bins_merged = TimeEnergyBins.merge_energy([bins2, self.bins])

    def test_merge_time(self):
        
        counts = [[10, 20, 30], [10, 20, 30]] 
        exposure = [1.0, 1.0]       
        # merge at high end
        tstart = [5.0, 6.0]
        tstop = [6.0, 7.0]
        bins2 = TimeEnergyBins(counts, tstart, tstop, exposure, 
                               self.bins.emin, self.bins.emax)
        bins_merged = TimeEnergyBins.merge_time([self.bins, bins2])
        self.assertListEqual(bins_merged.tstop.tolist(), [1.0, 2.0, 4.0, 5.0, 
                                                          6.0, 7.0])
        self.assertListEqual(bins_merged.tstart.tolist(), [0.0, 1.0, 3.0, 4.0,
                                                           5.0, 6.0])
        # flip order
        bins_merged = TimeEnergyBins.merge_time([bins2, self.bins])
        self.assertListEqual(bins_merged.tstop.tolist(), [1.0, 2.0, 4.0, 5.0, 
                                                          6.0, 7.0])
        self.assertListEqual(bins_merged.tstart.tolist(), [0.0, 1.0, 3.0, 4.0,
                                                           5.0, 6.0])
        # merge at low end
        tstart = [-2.0, -1.0]
        tstop = [-1.0, 0.0]
        bins2 = TimeEnergyBins(counts, tstart, tstop, exposure, 
                               self.bins.emin, self.bins.emax)
        bins_merged = TimeEnergyBins.merge_time([self.bins, bins2])
        self.assertListEqual(bins_merged.tstop.tolist(), [-1.0, 0.0, 1.0, 2.0, 
                                                          4.0, 5.0])
        self.assertListEqual(bins_merged.tstart.tolist(), [-2.0, -1.0, 0.0, 1.0,
                                                           3.0, 4.0])
                
        # overlapping merge (not allowed)
        tstart = [4.0, 5.0]
        tstop = [5.0, 6.0]
        bins2 = TimeEnergyBins(counts, tstart, tstop, exposure, 
                               self.bins.emin, self.bins.emax)
        with self.assertRaises(ValueError):
            bins_merged = TimeEnergyBins.merge_time([self.bins, bins2])
        with self.assertRaises(ValueError):
            bins_merged = TimeEnergyBins.merge_time([bins2, self.bins])

    def test_init_errors(self):
        with self.assertRaises(TypeError):
            TimeEnergyBins(0, self.bins.tstart, self.bins.tstop, 
                           self.bins.exposure, self.bins.emin, self.bins.emax)

        with self.assertRaises(TypeError):
            TimeEnergyBins(self.bins.counts[:,0], self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure, self.bins.emin,
                           self.bins.emax)
        
        with self.assertRaises(TypeError):
            TimeEnergyBins(self.bins.counts, self.bins.tstart[0], 
                           self.bins.tstop, self.bins.exposure, self.bins.emin,
                           self.bins.emax)

        with self.assertRaises(TypeError):
            TimeEnergyBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop[0], self.bins.exposure, 
                           self.bins.emin, self.bins.emax)

        with self.assertRaises(TypeError):
            TimeEnergyBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure[0], 
                           self.bins.emin, self.bins.emax)

        with self.assertRaises(TypeError):
            TimeEnergyBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure, 
                           self.bins.emin[0], self.bins.emax)

        with self.assertRaises(TypeError):
            TimeEnergyBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure, 
                           self.bins.emin, self.bins.emax[0])
        
        with self.assertRaises(ValueError):
            TimeEnergyBins(self.bins.counts, self.bins.tstart[1:], 
                           self.bins.tstop, self.bins.exposure, self.bins.emin,
                           self.bins.emax)

        with self.assertRaises(ValueError):
            TimeEnergyBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure[1:], 
                           self.bins.emin, self.bins.emax)

        with self.assertRaises(ValueError):
            TimeEnergyBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure[1:], 
                           self.bins.emin[1:], self.bins.emax)

        with self.assertRaises(ValueError):
            TimeEnergyBins(self.bins.counts.T, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure[1:], 
                           self.bins.emin, self.bins.emax)

        with self.assertRaises(TypeError):
            TimeEnergyBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure, 
                           self.bins.emin, self.bins.emax, quality=0)

        with self.assertRaises(ValueError):
            TimeEnergyBins(self.bins.counts, self.bins.tstart, 
                           self.bins.tstop, self.bins.exposure, 
                           self.bins.emin, self.bins.emax, quality=[0, 0])
         

class TestEventList(unittest.TestCase):
    
    def setUp(self):
        times = [0.706, 1.640, 3.185, 3.512, 5.540, 
                 7.790, 9.602, 9.726, 10.45, 10.61]
        chans = [0, 3, 2, 3, 2, 2, 3, 2, 2, 2]
        ebounds = Ebounds.from_bounds([10.0, 20.0, 40.0, 80.0], 
                                      [20.0, 40.0, 80.0, 160.0])
        self.ev = EventList(times=times, channels=chans, ebounds=ebounds)
        
    def test_channel_range(self):
        self.assertTupleEqual(self.ev.channel_range, (0, 3))
    
    def test_channels(self):
        self.assertListEqual(self.ev.channels.tolist(), 
                             [0, 3, 2, 3, 2, 2, 3, 2, 2, 2])

    def test_ebounds(self):
        self.assertIsInstance(self.ev.ebounds, Ebounds)
        
        ev = EventList(times=self.ev.times, channels=self.ev.channels)
        assert ev.ebounds is None
        
        with self.assertRaises(TypeError):
            ev.ebounds = [10.0, 20.0, 40.0, 80.0]
        
        ev.ebounds = self.ev.ebounds
        assert ev.ebounds == self.ev.ebounds
        
    def test_emax(self):
        self.assertEqual(self.ev.emax, 160.0)

    def test_emin(self):
        self.assertEqual(self.ev.emin, 10.0)

    def test_energy_range(self):
        self.assertTupleEqual(self.ev.energy_range, (10.0, 160.0))
    
    def test_num_chans(self):
        self.assertEqual(self.ev.num_chans, 4)
    
    def test_size(self):
        self.assertEqual(self.ev.size, 10)
    
    def test_time_range(self):
        self.assertTupleEqual(self.ev.time_range, (0.706, 10.61))

    def test_times(self):
        self.assertListEqual(self.ev.times.tolist(), 
                             [0.706, 1.640, 3.185, 3.512, 5.540, 
                              7.790, 9.602, 9.726, 10.45, 10.61])

    def test_bin(self):
        
        # full range, no deadtime, bin into 1 s starting at 0
        bins = self.ev.bin(bin_by_time, 1.0, time_ref=0.0, event_deadtime=0.0, 
                           overflow_deadtime=0.0)
        self.assertEqual(bins.num_chans, 4)
        self.assertEqual(bins.num_times, 11)
        counts = [[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 1, 1],
                  [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 1, 0],
                  [0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 2, 0]]
        for i in range(11):
             self.assertListEqual(bins.counts[i].tolist(), counts[i])
        self.assertListEqual(bins.exposure.tolist(), [1.0]*11)

        # full range, event deadtime, no overflow deadtime
        bins = self.ev.bin(bin_by_time, 1.0, time_ref=0.0, event_deadtime=0.1, 
                           overflow_deadtime=0.0)
        self.assertEqual(bins.num_chans, 4)
        self.assertEqual(bins.num_times, 11)
        self.assertListEqual(bins.exposure.tolist(), [0.9, 1.0, 1.0, 0.9, 1.0, 
                                                      0.9, 1.0, 0.9, 1.0, 0.9, 
                                                      0.8])        

        # full range, event deadtime, and overflow deadtime
        bins = self.ev.bin(bin_by_time, 1.0, time_ref=0.0, event_deadtime=0.1, 
                           overflow_deadtime=0.2)
        self.assertEqual(bins.num_chans, 4)
        self.assertEqual(bins.num_times, 11)
        self.assertListEqual(bins.exposure.tolist(), [0.9, 0.8, 1.0, 0.7, 1.0, 
                                                      0.9, 1.0, 0.9, 1.0, 0.7, 
                                                      0.8])

        # bin starting at tstart
        bins = self.ev.bin(bin_by_time, 1.0, time_ref=2.0, event_deadtime=0.0, 
                           overflow_deadtime=0.0, tstart=2.0)
        self.assertEqual(bins.num_chans, 4)
        self.assertEqual(bins.num_times, 9)
        counts = [[0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 0, 0], [0, 0, 1, 0], 
                  [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 1, 1], 
                  [0, 0, 2, 0]]
        for i in range(9):
             self.assertListEqual(bins.counts[i].tolist(), counts[i])
    
        # bin ending at tstop
        bins = self.ev.bin(bin_by_time, 1.0, time_ref=0.0, event_deadtime=0.0, 
                           overflow_deadtime=0.0, tstop=9.0)
        self.assertEqual(bins.num_chans, 4)
        self.assertEqual(bins.num_times, 9)
        counts = [[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 1, 1],
                  [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 1, 0],
                  [0, 0, 0, 0]]
        for i in range(9):
             self.assertListEqual(bins.counts[i].tolist(), counts[i])


        # bin starting at tstart and ending at tstop
        bins = self.ev.bin(bin_by_time, 1.0, time_ref=2.0, event_deadtime=0.0, 
                           overflow_deadtime=0.0, tstart=2.0, tstop=9.0)
        self.assertEqual(bins.num_chans, 4)
        self.assertEqual(bins.num_times, 7)
        counts = [[0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 0, 0], [0, 0, 1, 0], 
                  [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]]
        for i in range(7):
             self.assertListEqual(bins.counts[i].tolist(), counts[i])

        # bin outside of range
        bins = self.ev.bin(bin_by_time, 1.0, time_ref=20.0, event_deadtime=0.0, 
                           overflow_deadtime=0.0, tstart=20.0, tstop=22.0)
        self.assertEqual(bins.num_chans, 4)
        self.assertEqual(bins.num_times, 2)
        counts = [[0, 0, 0, 0], [0, 0, 0, 0]]
        for i in range(2):
             self.assertListEqual(bins.counts[i].tolist(), counts[i])
    
        # no ebounds
        ev = EventList(times=self.ev.times, channels=self.ev.channels)
        bins = ev.bin(bin_by_time, 1.0, time_ref=0.0, event_deadtime=0.0, 
                      overflow_deadtime=0.0)
        assert isinstance(bins, TimeChannelBins)
        self.assertEqual(bins.num_chans, 4)
        self.assertEqual(bins.num_times, 11)
        counts = [[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 1, 1],
                  [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 1, 0],
                  [0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 2, 0]]
        for i in range(11):
            self.assertListEqual(bins.counts[i].tolist(), counts[i])
        self.assertListEqual(bins.exposure.tolist(), [1.0]*11)
    
    def test_channel_slice(self):
        
        # slice at end of range
        ev2 = self.ev.channel_slice(2, 5)
        self.assertEqual(ev2.size, 9)        
        self.assertTupleEqual(ev2.channel_range, (2,3))        

        # slice at beginning of range
        ev2 = self.ev.channel_slice(0, 1)
        self.assertEqual(ev2.size, 1)        
        self.assertTupleEqual(ev2.channel_range, (0,0))        

        # slice single channel
        ev2 = self.ev.channel_slice(2, 2)
        self.assertEqual(ev2.size, 6)        
        self.assertTupleEqual(ev2.channel_range, (2,2))        

        # slice outside range
        ev2 = self.ev.channel_slice(5, 10)
        self.assertEqual(ev2.size, 0)        
        self.assertIsNone(ev2.channel_range)       

    def test_count_spectrum(self):
        # no deadtime
        bins = self.ev.count_spectrum()
        self.assertEqual(bins.size, 4)  
        self.assertListEqual(bins.counts.tolist(), [1, 0, 6, 3])    
        self.assertListEqual(bins.lo_edges.tolist(), [10.0, 20.0, 40.0, 80.0])
        self.assertListEqual(bins.hi_edges.tolist(), [20.0, 40.0, 80.0, 160.0])
        self.assertListEqual(bins.exposure.tolist(), [9.904]*4)
        
        # deadtime
        bins = self.ev.count_spectrum(event_deadtime=0.1, overflow_deadtime=0.2)
        self.assertEqual(bins.size, 4)  
        self.assertListEqual(bins.counts.tolist(), [1, 0, 6, 3])    
        self.assertListEqual(bins.exposure.tolist(), [8.604]*4)
        
        # no ebounds
        ev = EventList(times=self.ev.times, channels=self.ev.channels)
        bins = ev.count_spectrum()
        assert isinstance(bins, ChannelBins)
        self.assertEqual(bins.size, 4)  
        self.assertListEqual(bins.counts.tolist(), [1, 0, 6, 3])
        self.assertListEqual(bins.chan_nums.tolist(), [0, 1, 2, 3])  
        self.assertListEqual(bins.exposure.tolist(), [9.904]*4)
        
    def test_energy_slice(self):
        
        # slice at end of range
        ev2 = self.ev.energy_slice(50., 300.0)
        self.assertEqual(ev2.size, 9)        
        self.assertTupleEqual(ev2.channel_range, (2,3))        

        # slice at beginning of range
        ev2 = self.ev.energy_slice(15.0, 25.0)
        self.assertEqual(ev2.size, 1)        
        self.assertTupleEqual(ev2.channel_range, (0,0))        

        # slice single channel
        ev2 = self.ev.energy_slice(50.0, 50.0)
        self.assertEqual(ev2.size, 6)        
        self.assertTupleEqual(ev2.channel_range, (2,2))        

        # slice outside range
        ev2 = self.ev.channel_slice(300.0, 1000.0)
        self.assertEqual(ev2.size, 0)        
        self.assertIsNone(ev2.channel_range)       

    def test_get_exposure(self):
        
        # full range, no deadtime
        self.assertEqual(self.ev.get_exposure(), 9.904)
        
        # full range, event deadtime, no overflow deadtime
        self.assertEqual(self.ev.get_exposure(event_deadtime=0.1), 9.204)
        
        # full range, event and overflow deadtime
        self.assertEqual(self.ev.get_exposure(event_deadtime=0.1,
                                              overflow_deadtime=0.2), 8.604)
        
        # one time range
        self.assertEqual(self.ev.get_exposure(time_ranges=[1.0, 3.0]), 2.0)

        # two time ranges
        self.assertEqual(self.ev.get_exposure(time_ranges=[[1.0, 3.0],
                                                           [5.0, 7.0]]), 4.0)
    
    def test_rebin_energy(self):
        
        ev2 = self.ev.rebin_energy(combine_by_factor, 2)
        self.assertEqual(ev2.size, 10)
        self.assertEqual(ev2.num_chans, 2)
        self.assertListEqual(ev2.channels.tolist(), 
                             [0, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    
    def test_sort(self):
        self.ev.sort_channels()
        self.assertListEqual(self.ev.channels.tolist(), 
                             [0, 2, 2, 2, 2, 2, 2, 3, 3, 3])
        self.assertListEqual(self.ev.times.tolist(),
                            [0.706, 3.185, 5.540, 7.790, 9.726, 10.45, 10.61,
                             1.640, 3.512, 9.602])

        self.ev.sort_time()
        self.assertListEqual(self.ev.channels.tolist(), 
                             [0, 3, 2, 3, 2, 2, 3, 2, 2, 2])
        self.assertListEqual(self.ev.times.tolist(),
                            [0.706, 1.640, 3.185, 3.512, 5.540, 
                             7.790, 9.602, 9.726, 10.45, 10.61])

    def test_time_slice(self):
        
        # slice at end of range
        ev2 = self.ev.time_slice(9.0, 15.0)
        self.assertEqual(ev2.size, 4)        
        self.assertTupleEqual(ev2.time_range, (9.602, 10.61))        

        # slice at beginning of range
        ev2 = self.ev.time_slice(0.0, 4.0)
        self.assertEqual(ev2.size, 4)        
        self.assertTupleEqual(ev2.time_range, (0.706, 3.512))        

        # slice middle
        ev2 = self.ev.time_slice(5.0, 8.0)
        self.assertEqual(ev2.size, 2)        
        self.assertTupleEqual(ev2.time_range, (5.540, 7.790))        

        # slice outside range
        ev2 = self.ev.time_slice(15., 20.)
        self.assertEqual(ev2.size, 0)        
        self.assertIsNone(ev2.time_range)       

    def test_merge(self):
        times = [2.060, 2.971, 4.688, 4.938, 6.870, 
                 7.536, 8.200, 9.456, 10.91, 12.43]
        chans = [0, 0, 1, 3, 0, 1, 2, 2, 3, 0]
        ev2 = EventList(times=times, channels=chans, ebounds=self.ev.ebounds)

        # force unique
        ev_merged = EventList.merge([self.ev, ev2], sort=False)
        self.assertEqual(ev_merged.size, 20)
        self.assertListEqual(ev_merged.times.tolist(), 
                             [0.706, 1.640, 2.060, 2.971, 3.185, 3.512, 4.688, 
                              4.938, 5.540, 6.870, 7.536, 7.790, 8.200, 9.456, 
                              9.602, 9.726, 10.45, 10.61, 10.91, 12.43])
        
        # flip order
        ev_merged = EventList.merge([ev2, self.ev], sort=False)
        self.assertListEqual(ev_merged.times.tolist(), 
                             [0.706, 1.640, 2.060, 2.971, 3.185, 3.512, 4.688, 
                              4.938, 5.540, 6.870, 7.536, 7.790, 8.200, 9.456, 
                              9.602, 9.726, 10.45, 10.61, 10.91, 12.43])
        
        # do not force unique and sort
        ev_merged = EventList.merge([ev2, self.ev], sort=True, force_unique=False)
        self.assertEqual(ev_merged.size, 12)
        self.assertListEqual(ev_merged.times.tolist(), 
                             [0.706, 1.640, 3.185, 3.512, 5.540, 7.790, 9.602, 
                              9.726, 10.45, 10.61, 10.91, 12.43])
        
        # add a duplicate time
        times.append(3.185)
        chans.append(2)
        ev2 = EventList(times=times, channels=chans, ebounds=self.ev.ebounds)
        ev_merged = EventList.merge([self.ev, ev2], sort=False)
        self.assertEqual(ev_merged.size, 20)
        self.assertListEqual(ev_merged.times.tolist(), 
                             [0.706, 1.640, 2.060, 2.971, 3.185, 3.512, 4.688, 
                              4.938, 5.540, 6.870, 7.536, 7.790, 8.200, 9.456, 
                              9.602, 9.726, 10.45, 10.61, 10.91, 12.43])

    def test_init_errors(self):
        with self.assertRaises(ValueError):
            EventList(times=self.ev.times[1:], channels=self.ev.channels)
        
        with self.assertRaises(TypeError):
            EventList(times=self.ev.times, channels=self.ev.channels, ebounds=1)
        
        with self.assertRaises(AttributeError):
            ev = EventList(times=self.ev.time, channels=self.ev.channels)
            ev.energy_slice(50.0, 300.0)


class TestResponseMatrix(unittest.TestCase):
    
    def setUp(self):
        matrix = np.diag([0.1, 0.2, 0.3, 0.2, 0.1, 0.1])
        emin = [10., 20., 40., 80., 160., 320.]
        emax = [20., 40., 80., 160., 320., 640.]
        chanlo = [10.1, 20.1, 40.1, 80.1, 160.1, 320.1]
        chanhi = [20.1, 40.1, 80.1, 160.1, 320.1, 640.1]
        self.rsp = ResponseMatrix(matrix, emin, emax, chanlo, chanhi)
        
    def test_channel_centroids(self):
        vals = [14.25, 28.39, 56.67, 113.24, 226.38, 452.65]
        for i in range(4):
            self.assertAlmostEqual(self.rsp.channel_centroids[i], vals[i],
                                   places=2)
    
    def test_channel_widths(self):
        vals = [10.0, 20.0, 40.0, 80.0, 160., 320.]
        for i in range(4):
            self.assertAlmostEqual(self.rsp.channel_widths[i], vals[i],
                                   places=2)
    
    def test_ebounds(self):
        self.assertIsInstance(self.rsp.ebounds, Ebounds)
        self.assertTupleEqual(self.rsp.ebounds.range, (10.1, 640.1))
    
    def test_matrix(self):
        vals = [[0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.2, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.3, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.2, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.1, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.1]]
        for i in range(6):
            self.assertListEqual(self.rsp.matrix[i,:].tolist(), vals[i])

    def test_num_chans(self):
        self.assertEqual(self.rsp.num_chans, 6)
    
    def test_num_ebins(self):
        self.assertEqual(self.rsp.num_ebins, 6)

    def test_photon_bins(self):
        self.assertIsInstance(self.rsp.photon_bins, Ebounds)
        self.assertTupleEqual(self.rsp.photon_bins.range, (10., 640.0))

    def test_photon_bin_centroids(self):
        vals = [14.14, 28.28, 56.57, 113.14, 226.27, 452.55]
        for i in range(4):
            self.assertAlmostEqual(self.rsp.photon_bin_centroids[i], vals[i],
                                   places=2)
 
    def test_photon_bin_widths(self):
        vals = [10.0, 20.0, 40.0, 80.0, 160., 320.]
        for i in range(4):
            self.assertAlmostEqual(self.rsp.photon_bin_widths[i], vals[i],
                                   places=2)
   
    def test_channel_effective_area(self):
        effarea = self.rsp.channel_effective_area()
        self.assertListEqual(effarea.counts.tolist(), 
                             [0.1, 0.2, 0.3, 0.2, 0.1, 0.1])
    
    def test_effective_area(self):
        self.assertAlmostEqual(self.rsp.effective_area(10.0), 0.1)
        self.assertAlmostEqual(self.rsp.effective_area(40.0), 0.3)
        with self.assertRaises(ValueError):
            self.rsp.effective_area(1.0)
        with self.assertRaises(ValueError):
            self.rsp.effective_area(1000.0)
    
    def test_fold_spectrum(self):
        def zeroth_poly(c, x):
            return np.full(x.shape, c[0])
        
        counts = self.rsp.fold_spectrum(zeroth_poly, (1.0,))
        self.assertListEqual(counts.tolist(), [1., 4., 12., 16., 16., 32.])

        counts = self.rsp.fold_spectrum(zeroth_poly, (2.0,))
        self.assertListEqual(counts.tolist(), [2., 8., 24., 32., 32., 64.])
        
        # channel mask
        mask = np.array([False, True, True, True, True, False])
        counts = self.rsp.fold_spectrum(zeroth_poly, (1.0,), channel_mask=mask)
        self.assertListEqual(counts.tolist(), [4., 12., 16., 16.])

    def test_photon_effective_area(self):
        effarea = self.rsp.photon_effective_area()
        self.assertListEqual(effarea.counts.tolist(), 
                             [0.1, 0.2, 0.3, 0.2, 0.1, 0.1])
    
    def test_rebin(self):
        # rebin by factor
        rsp2 = self.rsp.rebin(factor=2)
        self.assertListEqual(rsp2.ebounds.low_edges(), [10.1, 40.1, 160.1])
        self.assertListEqual(rsp2.ebounds.high_edges(), [40.1, 160.1, 640.1])
        vals = [[0.1, 0.0, 0.0],
                [0.2, 0.0, 0.0],
                [0.0, 0.3, 0.0],
                [0.0, 0.2, 0.0],
                [0.0, 0.0, 0.1],
                [0.0, 0.0, 0.1]]
        for i in range(6):
            self.assertListEqual(rsp2.matrix[i,:].tolist(), vals[i])
        
        # rebin by edge indices
        rsp2 = self.rsp.rebin(edge_indices=[0,2,4,6])
        self.assertListEqual(rsp2.ebounds.low_edges(), [10.1, 40.1, 160.1])
        self.assertListEqual(rsp2.ebounds.high_edges(), [40.1, 160.1, 640.1])
        for i in range(6):
            self.assertListEqual(rsp2.matrix[i,:].tolist(), vals[i])

        # invalid input
        with self.assertRaises(ValueError):
            self.rsp.rebin()
        
        # invalid factor
        with self.assertRaises(ValueError):
            self.rsp.rebin(factor=0)

        # invalid factor
        with self.assertRaises(TypeError):
            self.rsp.rebin(factor='')

        # invalid factor
        with self.assertRaises(ValueError):
            self.rsp.rebin(factor=5)
            
        # invalid edge indices
        with self.assertRaises(ValueError):
            self.rsp.rebin(edge_indices=[-1, 5])

        # invalid edge indices
        with self.assertRaises(ValueError):
            self.rsp.rebin(edge_indices=[0, 10])

    def test_resample(self):
        rsp2 = self.rsp.resample(num_photon_bins=6)
        vals = [10., 20., 40., 80., 160., 320.]
        for i in range(6):
            self.assertAlmostEqual(rsp2.photon_bins.low_edges()[i], vals[i])
        vals = [20., 40., 80., 160., 320., 640.]
        for i in range(6):
            self.assertAlmostEqual(rsp2.photon_bins.high_edges()[i], vals[i])

        rsp2 = self.rsp.resample(num_photon_bins=3)
        vals = [10., 40., 160.]
        for i in range(3):
            self.assertAlmostEqual(rsp2.photon_bins.low_edges()[i], vals[i])
        vals = [40., 160., 640.]
        for i in range(3):
            self.assertAlmostEqual(rsp2.photon_bins.high_edges()[i], vals[i])
        
        rsp2 = self.rsp.resample(photon_bin_edges=[10.0, 40.0, 160.0, 640.])
        vals = [10., 40., 160.]
        for i in range(3):
            self.assertAlmostEqual(rsp2.photon_bins.low_edges()[i], vals[i])
        vals = [40., 160., 640.]
        for i in range(3):
            self.assertAlmostEqual(rsp2.photon_bins.high_edges()[i], vals[i])

        with self.assertRaises(ValueError):
            rsp2 = self.rsp.resample()        

        with self.assertRaises(TypeError):
            rsp2 = self.rsp.resample(num_photon_bins=-1)
        
        with self.assertRaises(TypeError):
            rsp2 = self.rsp.resample(num_photon_bins='')

        with self.assertRaises(ValueError):
            rsp2 = self.rsp.resample(photon_bin_edges=[0.0, 50.0])

        with self.assertRaises(ValueError):
            rsp2 = self.rsp.resample(photon_bin_edges=[50.0, 1000.0])

    def test_init_errors(self):
        matrix = np.diag([0.1, 0.2, 0.3, 0.2, 0.1, 0.1])
        emin = [10., 20., 40., 80., 160., 320.]
        emax = [20., 40., 80., 160., 320., 640.]
        chanlo = [10.1, 20.1, 40.1, 80.1, 160.1, 320.1]
        chanhi = [20.1, 40.1, 80.1, 160.1, 320.1, 640.1]

        with self.assertRaises(TypeError):
            ResponseMatrix(0, emin, emax, chanlo, chanhi)

        with self.assertRaises(TypeError):
            ResponseMatrix(matrix[:,0], emin, emax, chanlo, chanhi)

        with self.assertRaises(TypeError):
            ResponseMatrix(matrix, emin[0], emax, chanlo, chanhi)

        with self.assertRaises(TypeError):
            ResponseMatrix(matrix, emin, emax[0], chanlo, chanhi)

        with self.assertRaises(TypeError):
            ResponseMatrix(matrix, emin, emax, chanlo[0], chanhi)

        with self.assertRaises(TypeError):
            ResponseMatrix(matrix, emin, emax, chanlo, chanhi[0])

        with self.assertRaises(ValueError):
            ResponseMatrix(matrix, emin[1:], emax, chanlo, chanhi)

        with self.assertRaises(ValueError):
            ResponseMatrix(matrix, emin, emax, chanlo[1:], chanhi)

        with self.assertRaises(ValueError):
            ResponseMatrix(matrix[1:,:], emin, emax, chanlo[1:], chanhi)


class TestParameter(unittest.TestCase):
    
    def setUp(self):
        self.param = Parameter(100.0, (5.0, 10.0), name='Epeak',
                               units='keV', support=(0.0, np.inf))
    
    def test_name(self):
        self.assertEqual(self.param.name, 'Epeak')

    def test_support(self):
        self.assertTupleEqual(self.param.support, (0.0, np.inf))

    def test_uncertainty(self):
        self.assertTupleEqual(self.param.uncertainty, (5.0, 10.0))
        
        param = Parameter(100.0, 5.0)
        self.assertTupleEqual(param.uncertainty, (5.0, 5.0))

    def test_units(self):
        self.assertEqual(self.param.units, 'keV')        

    def test_value(self):
        self.assertEqual(self.param.value, 100.0)        

    def test_one_sigma_range(self):
        self.assertTupleEqual(self.param.one_sigma_range(), (95.0, 110.0))

    def test_to_fits_value(self):
        self.assertTupleEqual(self.param.to_fits_value(), (100.0, 10.0, 5.0))
        
        param = Parameter(100.0, 5.0)
        self.assertTupleEqual(param.to_fits_value(), (100.0, 5.0, 5.0))

    def test_valid_value(self):
        self.assertTrue(self.param.valid_value())
        
        param = Parameter(100.0, 5.0, support=(200.0, np.inf))
        self.assertFalse(param.valid_value())

    def test_init_errors(self):
        with self.assertRaises(ValueError):
            Parameter(100.0, (5.0, 10.0, 0.0))

        with self.assertRaises(TypeError):
            Parameter(100.0, '5.0')
        

if __name__ == '__main__':
    unittest.main()
