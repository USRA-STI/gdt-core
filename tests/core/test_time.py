#  CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
#  Contract No.: CA 80MSFC17M0022
#  Contractor Name: Universities Space Research Association
#  Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
#  Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
#  Developed by: William Cleveland and Adam Goldstein
#                Universities Space Research Association
#                Science and Technology Institute
#                https://sti.usra.edu
#
#  Developed by: Daniel Kocevski
#                National Aeronautics and Space Administration (NASA)
#                Marshall Space Flight Center
#                Astrophysics Branch (ST-12)
#
#   Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
#   in compliance with the License. You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software distributed under the License
#  is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
#  implied. See the License for the specific language governing permissions and limitations under the
#  License.
#
import unittest
import numpy as np
import astropy.units as u
from astropy.time import Time, TimeDelta
from gdt.core.time import time_range, truncate_hour, hours_range_from


class TestTimeRange(unittest.TestCase):

    def test_start_stop_default_step(self):
        start = Time('2022-09-09 12:00:00')
        stop = Time('2022-09-09 13:00:00')

        expected = Time([start + TimeDelta(x, format='sec') for x in range(3600)])

        result = time_range(start, stop)
        self.assertEqual(len(result), 3600)
        self.assertTrue(np.all(result.unix == expected.unix))

    def test_start_stop_5sec_step_primitive(self):
        start = Time('2022-09-09 12:00:00')
        stop = Time('2022-09-09 13:00:00')

        expected = Time([start + TimeDelta(x, format='sec') for x in range(0, 3600, 5)])

        result = time_range(start, stop, 5)
        self.assertEqual(len(result), 720)
        self.assertTrue(np.all(result.unix == expected.unix))

    def test_start_stop_15sec_step_time_delta(self):
        start = Time('2022-09-09 12:00:00')
        stop = Time('2022-09-09 13:00:00')

        expected = Time([start + TimeDelta(x, format='sec') for x in range(0, 3600, 15)])

        result = time_range(start, stop, TimeDelta(15, format='sec'))
        self.assertEqual(len(result), 240)
        self.assertTrue(np.all(result.unix == expected.unix))

    def test_start_stop_20min_step_quantity(self):
        start = Time('2022-09-09 12:00:00')
        stop = Time('2022-09-09 13:00:00')

        expected = Time([start + TimeDelta(x, format='sec') for x in range(0, 3600, 1200)])

        result = time_range(start, stop, 20 * u.min)
        self.assertEqual(len(result), 3)
        self.assertTrue(np.all(result.unix == expected.unix))

    def test_start_10min_interval_1min_step_primitive(self):
        start = Time('2022-09-09 12:00:00')

        expected = Time([start + TimeDelta(x, format='sec') for x in range(0, 600, 60)])

        result = time_range(start, 600, 1 * u.min)
        self.assertEqual(len(result), 10)
        self.assertTrue(np.all(result.unix == expected.unix))

    def test_start_15min_interval_1min_step_time_delta(self):
        start = Time('2022-09-09 12:00:00')

        expected = Time([start + TimeDelta(x, format='sec') for x in range(0, 900, 60)])

        result = time_range(start, TimeDelta(900, format='sec'), 1 * u.min)
        self.assertEqual(len(result), 15)
        self.assertTrue(np.all(result.unix == expected.unix))

    def test_start_6hr_interval_30min_step_quantity(self):
        start = Time('2022-09-09 12:00:00')

        expected = Time([start + TimeDelta(x, format='sec') for x in range(0, 21600, 1800)])

        result = time_range(start, 6 * u.hr, 30 * u.min)
        self.assertEqual(len(result), 12)
        self.assertTrue(np.all(result.unix == expected.unix))


class TestTruncateHour(unittest.TestCase):

    def test_truncate(self):
        dval = Time('2021-12-15 10:45:30', format='iso', scale='utc')
        self.assertEqual(truncate_hour(dval), Time('2021-12-15 10:00:00', format='iso', scale='utc'))


class TestHourRangeFrom(unittest.TestCase):

    def test_hours_from(self):
        dval = Time('2021-12-15 10:45:30', format='iso', scale='utc')
        expected =  Time([
            Time('2021-12-15 05:00:00', format='iso', scale='utc'),
            Time('2021-12-15 06:00:00', format='iso', scale='utc'),
            Time('2021-12-15 07:00:00', format='iso', scale='utc'),
            Time('2021-12-15 08:00:00', format='iso', scale='utc'),
            Time('2021-12-15 09:00:00', format='iso', scale='utc'),
        ])

        hours = hours_range_from(5, from_time=dval)
        for h, e in zip(hours, expected):
            self.assertEqual(h.to_datetime(), e.to_datetime())

