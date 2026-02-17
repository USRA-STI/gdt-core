# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Developed by: Joshua Wood
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Developed by: Jacob Smith
#               University of Alabama in Huntsville
#               Center for Space Plasma and Aeronomic Research
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
import numpy as np
import unittest

from astropy import units as u

from gdt.core.tte import PhotonList
from gdt.core.plot.lightcurve import Lightcurve
from gdt.core.data_primitives import Ebounds, Gti, EventList
from gdt.core.binning.unbinned import bin_by_time
from gdt.core.trigger import Trigger, TriggerAlgorithm, SlidingWindowMethod


class TestTrigger(unittest.TestCase):

    @classmethod 
    def setUpClass(cls):

        ebounds = Ebounds.from_bounds(
            [0, 1, 2, 3, 4, 5, 6, 7],
            [1, 2, 3, 4, 5, 6, 7, 8])

        # test TTE with flat background
        times = np.linspace(-50, 50, 10000)
        energy = np.arange(times.size) % 8
        cls.tte = PhotonList.from_data(EventList(times, energy, ebounds))

        # test TTE with 1 second transient
        times = np.concatenate([times, np.linspace(-0.4999, 0.4999, 200)])
        energy = np.arange(times.size) % 8
        cls.tte_transient = PhotonList.from_data(EventList(times, energy, ebounds))

        # test algorithms
        cls.algs = {
            1: TriggerAlgorithm(1000, 0, [3, 4], 4.5),
            2: TriggerAlgorithm(1000, 500, [3, 4], 4.5),
            3: TriggerAlgorithm(1000, 500, [2, 2], 6.0)
        }

        # initialize trigger method class
        cls.trigger_method = SlidingWindowMethod(
            algorithms=cls.algs, resolution=500,
            background_window=10000, background_offset=4000,
            channel_edges=[0, 1, 2, 3, 4, 5, 6, 7, 8],
            det_names=["det0", "det1", "det2"])
        cls.trigger_method.prepare_data([cls.tte_transient, cls.tte_transient, cls.tte])

        # create a set of triggers for use with plotting methods
        cls.triggers = cls.trigger_method.apply_trigger()

    def test_apply_trigger(self):

        triggers = self.trigger_method.apply_trigger()
        self.assertEqual(len(triggers), 4)
        self.assertEqual(triggers[0].alg_num, 1)
        self.assertEqual(triggers[0].alg.timescale, 1000)
        self.assertEqual(triggers[0].time, 0.0)
        self.assertEqual(triggers[0].triggered_det_names, ['det0', 'det1'])

        # check again with holdoff and debug options
        triggers_with_holdoff = self.trigger_method.apply_trigger(holdoff=300, debug=True)
        self.assertEqual(len(triggers_with_holdoff), 1)
        self.assertEqual(triggers[0].alg_num, triggers_with_holdoff[0].alg_num)
        self.assertEqual(triggers[0].alg.timescale, triggers_with_holdoff[0].alg.timescale)
        self.assertEqual(triggers[0].time, triggers_with_holdoff[0].time)
        self.assertEqual(triggers[0].triggered_det_names, triggers_with_holdoff[0].triggered_det_names)

    def test_lightcurve_plot(self):

        # single detector
        fig, axes, lcplots = self.trigger_method.lightcurve_plot(
            self.triggers[0], detectors=[0])
        # all detectors
        fig, axes, lcplots = self.trigger_method.lightcurve_plot(
            self.triggers[0])

    def test_waterfall_plot(self):

        rects = self.trigger_method.waterfall_plot(self.triggers)
        rects0 = self.trigger_method.waterfall_plot(
            [self.triggers[0]], facecolor='#FA8072', edgecolor='#B2342D', linewidth=3)

    def test_init_errors(self):

        # a poorly formated algorithms dictionary
        with self.assertRaises(ValueError):
            SlidingWindowMethod(
                algorithms=0, channel_edges=[0, 1], det_names=["det0"],
                resolution=500, background_window=10000, background_offset=1000)
        # algorithm timescale doesn't match time bin resolution
        with self.assertRaises(ValueError):
            algs = {1:TriggerAlgorithm(timescale=17, offset=0, channels=[0, 0], threshold=5.0)}
            SlidingWindowMethod(
                algorithms=algs, channel_edges=[0, 1], det_names=["det0"],
                resolution=500, background_window=10000, background_offset=1000)
        # algorithm offset doesn't match time bin resolution
        with self.assertRaises(ValueError):
            algs = {1:TriggerAlgorithm(timescale=1000, offset=17, channels=[0, 0], threshold=5.0)}
            SlidingWindowMethod(
                algorithms=algs, channel_edges=[0, 1], det_names=["det0"],
                resolution=500, background_window=10000, background_offset=1000)
        # background window doesn't match time bin resolution
        with self.assertRaises(ValueError):
            SlidingWindowMethod(
                algorithms=self.algs, channel_edges=[0, 1], det_names=["det0"],
                resolution=500, background_window=1, background_offset=1000)
        # background offset doesn't match time bin resolution
        with self.assertRaises(ValueError):
            SlidingWindowMethod(
                algorithms=self.algs, channel_edges=[0, 1], det_names=["det0"],
                resolution=500, background_window=10000, background_offset=1)

    def test_exec_errors(self):

        # applying trigger without preparing data
        with self.assertRaises(ValueError):
            trigger_method = SlidingWindowMethod(
                algorithms=self.algs, channel_edges=[0, 1], det_names=["det0"],
                resolution=500, background_window=10000, background_offset=1000)
            trigger_method.apply_trigger()
        # lightcurve where the requested time bin resolution is too small
        with self.assertRaises(ValueError):
            self.trigger_method.lightcurve_plot(self.triggers[0], resolution=1)
        # lightcurve without detectors to plot
        with self.assertRaises(ValueError):
            self.trigger_method.lightcurve_plot(self.triggers[0], detectors=[])
        # lightcurve without data preparation step
        with self.assertRaises(ValueError):
            trigger_method = SlidingWindowMethod(
                algorithms=self.algs, channel_edges=[0, 1], det_names=["det0"],
                resolution=500, background_window=10000, background_offset=1000)
            trigger_method.lightcurve_plot(self.triggers[0])

        # tte data without all detectors present
        with self.assertRaises(ValueError):
            self.trigger_method.prepare_data([self.tte])

    def test_apply_holdoff(self):

        # prepare some test triggers
        triggered_det = np.array([True, True])
        det_names = self.trigger_method.det_names
        alg = TriggerAlgorithm(timescale=16, offset=0, channels=[3,4], threshold=5.0)
        triggers = [
            Trigger(1, alg, 0.0, 5 * triggered_det.astype(float), triggered_det, det_names),
            # add another trigger with higher significance
            Trigger(1, alg, 0.0, 6 * triggered_det.astype(float), triggered_det, det_names),
            # add second trigger outside holdoff window
            Trigger(1, alg, 350.0, 6 * triggered_det.astype(float), triggered_det, det_names)
        ]

        selected = self.trigger_method.apply_holdoff(triggers, holdoff=300)
        self.assertEqual(len(selected), 2)
        self.assertEqual(selected[0].sig[0], 6.0)


class TestTriggerAlgorithm(unittest.TestCase):

    def test_alg(self):

        alg = TriggerAlgorithm(timescale=16, offset=16, channels=[3,4], threshold=5.0)
        self.assertEqual(alg.timescale_in_unit(u.second).value, 0.016)
        self.assertEqual(alg.offset_in_unit(u.second).value, 0.016)
        self.assertEqual(str(alg), "TriggerAlgorithm(timescale   16 ms, offset   16 ms, channels [3, 4], threshold 5.00 sigma)")
        self.assertEqual(alg.timescale, 16)
        self.assertEqual(alg.offset, 16)
        self.assertEqual(alg.channels, [3,4])
        self.assertEqual(alg.threshold, 5.0)

    def test_errors(self):

        # floating point timescale
        with self.assertRaises(ValueError):
            alg = TriggerAlgorithm(timescale=0.016, offset=16, channels=[3,4], threshold=5.0)
        # floating point offset
        with self.assertRaises(ValueError):
            alg = TriggerAlgorithm(timescale=16, offset=0.016, channels=[3,4], threshold=5.0)
        # lowest energy channel exceeds highest energy channel
        with self.assertRaises(ValueError):
            alg = TriggerAlgorithm(timescale=16, offset=16, channels=[4,3], threshold=5.0)
