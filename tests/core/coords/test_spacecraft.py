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
import astropy.coordinates.representation as r
import pytest
from astropy.coordinates import SkyCoord, GCRS, get_moon
from astropy.time import Time
from gdt.core.coords import SpacecraftFrame, SpacecraftAxes, Quaternion
from gdt.core.coords.spacecraft.axes import SpacecraftAxesAttribute

from astropy.coordinates import FunctionTransform, ICRS, frame_transform_graph
from gdt.core.coords.spacecraft.frame import icrs_to_spacecraft


class TestSpacecraftAxes(unittest.TestCase):

    def setUp(self):
        x_axis = SkyCoord(19.09, 8.04, frame='icrs', unit='deg')
        z_axis = SkyCoord(108.16, -6.54, frame='icrs', unit='deg')
        self.axes = SpacecraftAxes(x_pointing=x_axis, z_pointing=z_axis)

    def test_has_x(self):
        self.assertTrue(self.axes.has_x)

    def test_has_y(self):
        self.assertFalse(self.axes.has_y)

    def test_has_z(self):
        self.assertTrue(self.axes.has_z)

    def test_x_vector(self):
        self.assertListEqual(self.axes.x_vector.tolist(), [1.0, 0.0, 0.0])

    def test_y_vector(self):
        self.assertListEqual(self.axes.y_vector.tolist(), [0.0, 1.0, 0.0])

    def test_z_vector(self):
        self.assertListEqual(self.axes.z_vector.tolist(), [0.0, 0.0, 1.0])

    def test_pointing_vector(self):
        vals = [0.936, 0.324, 0.140]
        vec = self.axes.pointing_vector('x')
        for i in range(3):
            self.assertAlmostEqual(vec[i], vals[i], places=3)

        self.assertIsNone(self.axes.pointing_vector('y'))

        vals = [-0.310, 0.944, -0.114]
        vec = self.axes.pointing_vector('z')
        for i in range(3):
            self.assertAlmostEqual(vec[i], vals[i], places=3)

        axes = SpacecraftAxes(y_pointing=SkyCoord(19.09, 8.04, unit='deg'))
        vals = [0.936, 0.324, 0.140]
        vec = axes.pointing_vector('y')
        for i in range(3):
            self.assertAlmostEqual(vec[i], vals[i], places=3)

        with self.assertRaises(ValueError):
            self.axes.pointing_vector('hello')

    def test_init_errors(self):
        with self.assertRaises(TypeError):
            SpacecraftAxes(x_pointing=(0.0, 0.0))
        with self.assertRaises(TypeError):
            SpacecraftAxes(y_pointing=(0.0, 0.0))
        with self.assertRaises(TypeError):
            SpacecraftAxes(z_pointing=(0.0, 0.0))


class TestSpacecraftFrame(unittest.TestCase):

    def setUp(self):
        # Bare minimum frame
        quat = Quaternion.from_xyz_w(xyz=[0.0, 1.0, 0.0], w=1.0)
        self.min_frame = SpacecraftFrame(quaternion=quat)

        # Frame representing spacecraft position
        eic = r.CartesianRepresentation(-6320675.5, -1513143.1, 2313154.5, unit='m')
        t = Time('2016-12-31 23:58:56.740', format='iso', scale='utc')
        self.sc_frame = SpacecraftFrame(obsgeoloc=eic, obstime=t)

    def test_to_spacecraft_frame(self):
        coord = SkyCoord(0.0, 90.0, unit='deg', frame='icrs')
        transformed = coord.transform_to(self.min_frame)
        self.assertAlmostEqual(transformed.az.value[0], 180.0, places=5)
        self.assertAlmostEqual(transformed.el.value[0], 0.0, places=5)

    def test_to_equatorial(self):
        coord = SkyCoord(180.0, 0.0, unit='deg', frame=self.min_frame)
        self.assertAlmostEqual(coord.icrs.ra.value[0], 90.0, places=5)
        self.assertAlmostEqual(coord.icrs.dec.value[0], 90.0, places=5)

    def test_init_errors(self):
        with self.assertRaises(TypeError):
            SpacecraftFrame(quaternion=[0.0, 1.0, 0.0, 'a'])

    def test_height(self):
        self.assertAlmostEqual(self.sc_frame.earth_location.height.value, 522895.0, delta=0.5)

    def test_vector(self):
        self.assertListEqual(self.sc_frame.obsgeoloc.xyz.value.tolist(), [-6320675.5, -1513143.1, 2313154.5])

    def test_sc_gcrs(self):
        gcrs = self.sc_frame.sc_gcrs
        self.assertIsInstance(gcrs, GCRS)
        self.assertAlmostEqual(gcrs.ra.value, 193.46, places=2)
        self.assertAlmostEqual(gcrs.dec.value, 19.59, places=2)

    def test_geocenter(self):
        geocenter = self.sc_frame.geocenter
        self.assertAlmostEqual(geocenter.ra.value, 13.46, places=2)
        self.assertAlmostEqual(geocenter.dec.value, -19.59, places=2)

    def test_latitude(self):
        self.assertAlmostEqual(self.sc_frame.earth_location.lat.deg, 19.6, delta=0.3)

    def test_longitude(self):
        self.assertAlmostEqual(self.sc_frame.earth_location.lon.deg, 92.9, delta=0.3)

    def test_size(self):
        self.assertEqual(self.sc_frame.obsgeoloc.size, 1)

    def test_earth_angular_radius(self):
        self.assertAlmostEqual(self.sc_frame.earth_angular_radius.value, 67.55, places=2)

    def test_interpolate(self):
        eic = r.CartesianRepresentation(
            x=[-6330912.5, -6330020.5, -6329120.5, -6328213., -6327297.5,
               -6326374.5, -6325444., -6324506., -6323560., -6322606.,
               -6321644.5, -6320675.5, -6319699., -6318714., -6317721.5,
               -6316721.],
            y=[-1433237.5, -1440510.6, -1447781.9, -1455051.5, -1462319.1,
               -1469585.2, -1476849.5, -1484111.9, -1491372.5, -1498631.4,
               -1505888.2, -1513143.1, -1520396.5, -1527647.9, -1534897.2,
               -1542144.8],
            z=[2335938.2, 2333881.2, 2331821.5, 2329758.8, 2327692.8, 2325624.5,
               2323553.2, 2321479.2, 2319402.2, 2317322.2, 2315239.8, 2313154.5,
               2311066.2, 2308975., 2306881., 2304784.],
            unit='m')

        vel = r.CartesianRepresentation(
            x=[888.4821, 896.17035, 903.85846, 911.54565, 919.2348,
               926.91956, 934.604, 942.28564, 949.96826, 957.6493,
               965.3296, 973.00946, 980.68604, 988.36383, 996.0431,
               1003.7192],
            y=[-7273.888, -7272.1436, -7270.389, -7268.627, -7266.854,
               -7265.074, -7263.2847, -7261.487, -7259.6797, -7257.864,
               -7256.039, -7254.204, -7252.362, -7250.5103, -7248.649,
               -7246.779],
            z=[-2055.7422, -2058.5857, -2061.4268, -2064.2654, -2067.1033,
               -2069.9368, -2072.7678, -2075.596, -2078.422, -2081.2466,
               -2084.0679, -2086.8857, -2089.7017, -2092.5164, -2095.3289,
               -2098.1384],
            unit='m/s')

        quat = Quaternion([[-0.21510095, 0.01224809, 0.65282628, -0.7262227],
                           [-0.21533829, 0.01199657, 0.65275072, -0.72622449],
                           [-0.21557656, 0.01174484, 0.65267341, -0.72622739],
                           [-0.21581467, 0.01149171, 0.65259787, -0.72622861],
                           [-0.21605267, 0.01124013, 0.65252022, -0.72623156],
                           [-0.21629039, 0.01098749, 0.65244523, -0.72623203],
                           [-0.21652891, 0.01073432, 0.65236743, -0.72623463],
                           [-0.21676588, 0.01048162, 0.65229169, -0.72623567],
                           [-0.21700332, 0.01022866, 0.65221568, -0.72623663],
                           [-0.2172409, 0.00997442, 0.65213896, -0.72623803],
                           [-0.21747739, 0.00972049, 0.65206381, -0.72623817],
                           [-0.21771399, 0.00946618, 0.65198662, -0.72623994],
                           [-0.21795105, 0.00920997, 0.65191157, -0.72623951],
                           [-0.21818705, 0.00895516, 0.65183456, -0.72624095],
                           [-0.21842248, 0.00869974, 0.65175951, -0.72624065],
                           [-0.21865889, 0.0084436, 0.6516824, -0.72624172]])

        times = Time(['2016-12-31 23:58:45.740', '2016-12-31 23:58:46.740',
                      '2016-12-31 23:58:47.740', '2016-12-31 23:58:48.740',
                      '2016-12-31 23:58:49.740', '2016-12-31 23:58:50.740',
                      '2016-12-31 23:58:51.740', '2016-12-31 23:58:52.740',
                      '2016-12-31 23:58:53.740', '2016-12-31 23:58:54.740',
                      '2016-12-31 23:58:55.740', '2016-12-31 23:58:56.740',
                      '2016-12-31 23:58:57.740', '2016-12-31 23:58:58.740',
                      '2016-12-31 23:58:59.740', '2016-12-31 23:59:00.740'], format='iso', scale='utc')

        inter_times = Time(['2016-12-31 23:58:56.740', '2016-12-31 23:58:58.000',
                            '2016-12-31 23:58:59.000'], format='iso', scale='utc')
        sc_frame = SpacecraftFrame(obsgeoloc=eic, obsgeovel=vel, quaternion=quat, obstime=times)

        expected_lat = np.asarray([19.61267645, 19.58940036, 19.57090224])
        expected_lon = np.asarray([93.10018602, 93.1757562, 93.23571604])
        expected_height = np.asarray([522894.88816699, 522888.73016929, 522883.86069484])
        scpos_interp = sc_frame.at(inter_times)

        for i in range(len(inter_times)):
            self.assertAlmostEqual(scpos_interp[i].earth_location.lat.value, expected_lat[i], places=1)
            self.assertAlmostEqual(scpos_interp[i].earth_location.lon.value, expected_lon[i], places=1)
            self.assertAlmostEqual(scpos_interp[i].earth_location.height.value, expected_height[i], places=1)

    def test_interpolation_no_time(self):
        eic = r.CartesianRepresentation(
            x=[-6330912.5, -6330020.5, -6329120.5, -6328213., -6327297.5,
               -6326374.5, -6325444., -6324506., -6323560., -6322606.,
               -6321644.5, -6320675.5, -6319699., -6318714., -6317721.5,
               -6316721.],
            y=[-1433237.5, -1440510.6, -1447781.9, -1455051.5, -1462319.1,
               -1469585.2, -1476849.5, -1484111.9, -1491372.5, -1498631.4,
               -1505888.2, -1513143.1, -1520396.5, -1527647.9, -1534897.2,
               -1542144.8],
            z=[2335938.2, 2333881.2, 2331821.5, 2329758.8, 2327692.8, 2325624.5,
               2323553.2, 2321479.2, 2319402.2, 2317322.2, 2315239.8, 2313154.5,
               2311066.2, 2308975., 2306881., 2304784.],
            unit='m')

        sc_frame = SpacecraftFrame(obsgeoloc=eic)

        with pytest.raises(ValueError):
            sc_frame.at(Time.now())

    def test_interpolate_no_quat_no_vel(self):
        eic = r.CartesianRepresentation(
            x=[-6330912.5, -6330020.5, -6329120.5, -6328213., -6327297.5,
               -6326374.5, -6325444., -6324506., -6323560., -6322606.,
               -6321644.5, -6320675.5, -6319699., -6318714., -6317721.5,
               -6316721.],
            y=[-1433237.5, -1440510.6, -1447781.9, -1455051.5, -1462319.1,
               -1469585.2, -1476849.5, -1484111.9, -1491372.5, -1498631.4,
               -1505888.2, -1513143.1, -1520396.5, -1527647.9, -1534897.2,
               -1542144.8],
            z=[2335938.2, 2333881.2, 2331821.5, 2329758.8, 2327692.8, 2325624.5,
               2323553.2, 2321479.2, 2319402.2, 2317322.2, 2315239.8, 2313154.5,
               2311066.2, 2308975., 2306881., 2304784.],
            unit='m')

        times = Time(['2016-12-31 23:58:45.740', '2016-12-31 23:58:46.740',
                      '2016-12-31 23:58:47.740', '2016-12-31 23:58:48.740',
                      '2016-12-31 23:58:49.740', '2016-12-31 23:58:50.740',
                      '2016-12-31 23:58:51.740', '2016-12-31 23:58:52.740',
                      '2016-12-31 23:58:53.740', '2016-12-31 23:58:54.740',
                      '2016-12-31 23:58:55.740', '2016-12-31 23:58:56.740',
                      '2016-12-31 23:58:57.740', '2016-12-31 23:58:58.740',
                      '2016-12-31 23:58:59.740', '2016-12-31 23:59:00.740'], format='iso', scale='utc')

        inter_times = Time(['2016-12-31 23:58:56.740', '2016-12-31 23:58:58.000',
                            '2016-12-31 23:58:59.000'], format='iso', scale='utc')
        sc_frame = SpacecraftFrame(obsgeoloc=eic, obstime=times)

        expected_lat = np.asarray([19.61267645, 19.58940036, 19.57090224])
        expected_lon = np.asarray([93.10018602, 93.1757562, 93.23571604])
        expected_height = np.asarray([522894.88816699, 522888.73016929, 522883.86069484])
        scpos_interp = sc_frame.at(inter_times)

        for i in range(len(inter_times)):
            self.assertAlmostEqual(scpos_interp[i].earth_location.lat.value, expected_lat[i], places=1)
            self.assertAlmostEqual(scpos_interp[i].earth_location.lon.value, expected_lon[i], places=1)
            self.assertAlmostEqual(scpos_interp[i].earth_location.height.value, expected_height[i], places=1)

    def test_visibility(self):

        # Sun not visible, Moon is visible.
        time_1 = Time('2016-12-31 23:19:10.940', format='iso', scale='utc')
        sc_1 = SpacecraftFrame(obsgeoloc=r.CartesianRepresentation(5072027.5, 4579065.5, -1089665.6, unit='m'),
                               obstime=time_1)
        # Sun is visible, Moon is not visible.
        time_2 = Time('2016-12-31 23:55:25.740', format='iso', scale='utc')
        sc_2 = SpacecraftFrame(obsgeoloc=r.CartesianRepresentation(-6353991.5, 44461.45, 2687083.5, unit='m'),
                               obstime=time_2)

        self.assertTrue(sc_1.location_visible(get_moon(time_1)))
        self.assertFalse(sc_1.sun_visible)
        self.assertFalse(sc_2.location_visible(get_moon(time_2)))
        self.assertTrue(sc_2.sun_visible)

    def test_interp_velocity(self):
        eic = r.CartesianRepresentation(
            x=[-6330912.5, -6330020.5, -6329120.5, -6328213., -6327297.5,
               -6326374.5, -6325444., -6324506., -6323560., -6322606.,
               -6321644.5, -6320675.5, -6319699., -6318714., -6317721.5,
               -6316721.],
            y=[-1433237.5, -1440510.6, -1447781.9, -1455051.5, -1462319.1,
               -1469585.2, -1476849.5, -1484111.9, -1491372.5, -1498631.4,
               -1505888.2, -1513143.1, -1520396.5, -1527647.9, -1534897.2,
               -1542144.8],
            z=[2335938.2, 2333881.2, 2331821.5, 2329758.8, 2327692.8, 2325624.5,
               2323553.2, 2321479.2, 2319402.2, 2317322.2, 2315239.8, 2313154.5,
               2311066.2, 2308975., 2306881., 2304784.],
            unit='m')

        vel = r.CartesianRepresentation(
            x=[888.4821, 896.17035, 903.85846, 911.54565, 919.2348,
               926.91956, 934.604, 942.28564, 949.96826, 957.6493,
               965.3296, 973.00946, 980.68604, 988.36383, 996.0431,
               1003.7192],
            y=[-7273.888, -7272.1436, -7270.389, -7268.627, -7266.854,
               -7265.074, -7263.2847, -7261.487, -7259.6797, -7257.864,
               -7256.039, -7254.204, -7252.362, -7250.5103, -7248.649,
               -7246.779],
            z=[-2055.7422, -2058.5857, -2061.4268, -2064.2654, -2067.1033,
               -2069.9368, -2072.7678, -2075.596, -2078.422, -2081.2466,
               -2084.0679, -2086.8857, -2089.7017, -2092.5164, -2095.3289,
               -2098.1384],
            unit='m/s')

        quat = Quaternion([[-0.21510095, 0.01224809, 0.65282628, -0.7262227],
                           [-0.21533829, 0.01199657, 0.65275072, -0.72622449],
                           [-0.21557656, 0.01174484, 0.65267341, -0.72622739],
                           [-0.21581467, 0.01149171, 0.65259787, -0.72622861],
                           [-0.21605267, 0.01124013, 0.65252022, -0.72623156],
                           [-0.21629039, 0.01098749, 0.65244523, -0.72623203],
                           [-0.21652891, 0.01073432, 0.65236743, -0.72623463],
                           [-0.21676588, 0.01048162, 0.65229169, -0.72623567],
                           [-0.21700332, 0.01022866, 0.65221568, -0.72623663],
                           [-0.2172409, 0.00997442, 0.65213896, -0.72623803],
                           [-0.21747739, 0.00972049, 0.65206381, -0.72623817],
                           [-0.21771399, 0.00946618, 0.65198662, -0.72623994],
                           [-0.21795105, 0.00920997, 0.65191157, -0.72623951],
                           [-0.21818705, 0.00895516, 0.65183456, -0.72624095],
                           [-0.21842248, 0.00869974, 0.65175951, -0.72624065],
                           [-0.21865889, 0.0084436, 0.6516824, -0.72624172]])

        times = Time(['2016-12-31 23:58:45.740', '2016-12-31 23:58:46.740',
                      '2016-12-31 23:58:47.740', '2016-12-31 23:58:48.740',
                      '2016-12-31 23:58:49.740', '2016-12-31 23:58:50.740',
                      '2016-12-31 23:58:51.740', '2016-12-31 23:58:52.740',
                      '2016-12-31 23:58:53.740', '2016-12-31 23:58:54.740',
                      '2016-12-31 23:58:55.740', '2016-12-31 23:58:56.740',
                      '2016-12-31 23:58:57.740', '2016-12-31 23:58:58.740',
                      '2016-12-31 23:58:59.740', '2016-12-31 23:59:00.740'], format='iso', scale='utc')

        sc_frame = SpacecraftFrame(obsgeoloc=eic, obsgeovel=vel, quaternion=quat, obstime=times)

        # multiple times
        vel_times = times[1:-1]
        est_vel = sc_frame.estimate_velocity(vel_times, .1)
        orig_vel = vel[1:-1]
        diff = (np.abs(orig_vel.xyz.value - est_vel.xyz.value) / orig_vel.xyz.value) * 100

        # Verify the difference between telemetry and estimated velocity is less than 0.05 percent
        self.assertTrue(np.all(diff < 0.05))

        # single times
        vel_time = times[4]
        est_vel = sc_frame.estimate_velocity(vel_time, .1)
        orig_vel = vel[4]
        diff = (np.abs(orig_vel.xyz.value - est_vel.xyz.value) / orig_vel.xyz.value) * 100

        # Verify the difference between telemetry and estimated velocity is less than 0.05 percent
        self.assertTrue(np.all(diff < 0.05))

        # Test invalid time value
        with self.assertRaises(TypeError):
            sc_frame.estimate_velocity(12345678.0, 0.1)


class TestSpacecraftAxesAttribute(unittest.TestCase):

    def test_attribute_is_none(self):
        attr = SpacecraftAxesAttribute()

        value, converted = attr.convert_input(None)
        self.assertIsNone(value)
        self.assertFalse(converted)

    def test_attribute_is_bad_type(self):
        attr = SpacecraftAxesAttribute()

        with self.assertRaises(TypeError):
            value, converted = attr.convert_input(Time.now())

    def test_attribute_spacecraft_axes(self):
        attr = SpacecraftAxesAttribute()

        x_axis = SkyCoord(19.09, 8.04, frame='icrs', unit='deg')
        z_axis = SkyCoord(108.16, -6.54, frame='icrs', unit='deg')
        axes = SpacecraftAxes(x_pointing=x_axis, z_pointing=z_axis)

        value, converted = attr.convert_input(axes)
        self.assertTrue(isinstance(value, SpacecraftAxes))
        self.assertFalse(converted)
        self.assertEqual(axes, value)


class TestFrameTransfer(unittest.TestCase):

    def setUp(self):
        class TestFrame(SpacecraftFrame):
            pass
        quaternion_test = [-0.5636783622078214, 0.35243503385029756,
                            -0.34740741282879534, 0.6613352708008627]
        
        self.sc_frame = TestFrame(quaternion=quaternion_test)
        self.z_point = SkyCoord(30.271184, 6.667796, unit='deg')

        @frame_transform_graph.transform(FunctionTransform, ICRS, TestFrame)
        def icrs_to_test_frame(icrs_frame, test_frame):
            return icrs_to_spacecraft(icrs_frame, test_frame)

    def test_icrs_to_spacecraft(self):
        zaxis_test =self.z_point.transform_to(self.sc_frame)
        self.assertAlmostEqual(zaxis_test.el.value[0], 90.0, places=3)
       



