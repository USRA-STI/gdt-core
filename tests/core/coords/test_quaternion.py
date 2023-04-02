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
import pytest
import numpy as np
from astropy.time import Time
from scipy.spatial.transform import Rotation

from gdt.core.coords.quaternion import Quaternion, QuaternionAttribute


# Test scalar version of quaternions

@pytest.mark.parametrize("q1, scalar_first, expected", [[(1, 2, 3, 4), False, (1, 2, 3, 4)],
                                                        [(5, 6, 7, 8), True, (6, 7, 8, 5)]])
def test_to_xyz_w(q1, scalar_first, expected):
    quat = Quaternion(q1, scalar_first=scalar_first)
    assert quat == Quaternion(expected)
    assert np.all(quat.xyz == expected[:3])
    assert quat.w == expected[3]


@pytest.mark.parametrize("vector, scalar, expected, is_scalar", [
    [(1, 2, 3), 4, (1, 2, 3, 4), True],
    [(5, 6, 7), 8, (5, 6, 7, 8), True],
    [[(1, 2, 3), (5, 6, 7)], [4, 8], [(1, 2, 3, 4), (5, 6, 7, 8)], False]])
def test_from_xyz_w(vector, scalar, expected, is_scalar):
    quat = Quaternion.from_xyz_w(xyz=vector, w=scalar)
    assert quat == Quaternion(expected)
    assert quat.isscalar == is_scalar


@pytest.mark.parametrize("scalar_last, x, y, z, w", [[(1, 2, 3, 4), 1, 2, 3, 4]])
def test_component_properties(scalar_last, x, y, z, w):
    q = Quaternion(scalar_last)
    assert q.x == x
    assert q.y == y
    assert q.z == z
    assert q.w == w


@pytest.mark.parametrize("q, n", [[(1, 2, 3, 4), (-1, -2, -3, -4)],
                                  [(5, -6, 7, -8), (-5, 6, -7, 8)]])
def test_equals_and_negative(q, n):
    q_pos = Quaternion(q)
    q_neg = Quaternion(n)
    assert q_pos == q_pos
    assert q_pos != q_neg
    assert -q_pos == q_neg


@pytest.mark.parametrize("scalar_last, scalar_first", [[(1, 2, 3, 4), (4, 1, 2, 3)],
                                                       [(5, 6, 7, 8), (8, 5, 6, 7)],
                                                       [(9, 10, 11, 12), (12, 9, 10, 11)]])
def test_scalar_last_and_first(scalar_last, scalar_first):
    q = Quaternion(scalar_last)
    assert np.all(q.scalar_first == scalar_first)
    assert np.all(q.scalar_last == scalar_last)

    q = Quaternion(scalar_first, scalar_first=True)
    assert np.all(q.scalar_first == scalar_first)
    assert np.all(q.scalar_last == scalar_last)


@pytest.mark.parametrize("v1, v2, expected, is_scalar", [
    [(1, 2, 3), (4, 5, 6), (-3, 6, -3, 64.83291), True],
    [(7, 8, 9), (10, 11, 12), (-3, 6, -3, 532.10148), True],
    [[(1, 2, 3), (7, 8, 9)], [(4, 5, 6), (10, 11, 12)], [(-3, 6, -3, 64.83291), (-3, 6, -3, 532.10148)], False]])
def test_from_vectors(v1, v2, expected, is_scalar):
    q = Quaternion.from_vectors(v1, v2)
    assert q.round(5) == Quaternion(expected)
    assert q.isscalar == is_scalar


@pytest.mark.parametrize("q, expected", [[(9, 6, -4, 2), [[0.240876, 0.905109, -0.350365],
                                                          [0.671533, -0.416058, -0.613139],
                                                          [-0.70073, -0.0875912, -0.708029]]],
                                         [(-5, 8, -3, 0), [[-0.489796, -0.816327, 0.306122],
                                                           [-0.816327, 0.306122, -0.489796],
                                                           [0.306122, -0.489796, -0.816327]]],
                                         [(8, 13, 12, -16), [[0.0110585, 0.935229, -0.35387],
                                                             [-0.278041, 0.342812, 0.897314],
                                                             [0.960506, 0.0884676, 0.263823]]]])
def test_to_and_from_rotation(q, expected):
    quat = Quaternion(q)
    rot = quat.rotation.as_matrix()
    assert np.all(rot.round(decimals=6) == np.asarray(expected).round(decimals=6))

    rot = Rotation.from_matrix(expected)
    assert Quaternion.from_rotation(rot).equal_rotation(Quaternion(q).unit)


@pytest.mark.parametrize("q, expected", [[(9, 6, -4, 2), (-0.066, -0.044, 0.029, 0.015)],
                                         [(-5, 8, -3, 0), (0.051, -0.082, 0.031, 0.000)],
                                         [(8, 13, 12, -16), (-0.013, -0.021, -0.019, -0.025)]])
def test_inverse(q, expected):
    quat = Quaternion(q)
    assert quat.inverse.round(3) == Quaternion(expected)


@pytest.mark.parametrize("q, expected", [[(9, 6, -4, 2), 11.705],
                                         [(-5, 8, -3, 0), 9.899],
                                         [(8, 13, 12, -16), 25.159]])
def test_norm(q, expected):
    quat = Quaternion(q)
    assert round(quat.norm, 3) == round(expected, 3)


@pytest.mark.parametrize("q, expected", [[(9, 6, -4, 2), (0.769, 0.513, -0.342, 0.171)],
                                         [(-5, 8, -3, 0), (-0.505, 0.808, -0.303, 0)],
                                         [(8, 13, 12, -16), (0.318, 0.517, 0.477, -0.636)]])
def test_unit(q, expected):
    quat = Quaternion(q)
    assert quat.unit.round(3) == Quaternion(expected)


@pytest.mark.parametrize("q", [[(0.769, 0.513, -0.342, 0.171)],
                               [(-0.505, 0.808, -0.303, 0)],
                               [(0.318, 0.517, 0.477, -0.636)]])
def test_equal_rotation(q):
    q_pos = Quaternion(q)
    q_neg = -q_pos
    assert q_pos != q_neg
    assert q_pos.equal_rotation(q_neg)


@pytest.mark.parametrize("q1, q2, expected",
                         [[(-1, 3, 1, 2), (-4, 0, 1, 5), (-10, 12, 19, 5)],
                          [(-4, 0, 1, 5), (-1, 3, 1, 2), (-16, 18, -5, 5)],
                          [(8, 13, 12, -16), (-4, 10, 8, 15), (168, -77, 184, -434)],
                          [(-4, 10, 8, 15), (8, 13, 12, -16), (200, 147, -80, -434)]])
def test_multiply_quaternion(q1, q2, expected):
    quat1 = Quaternion(q1)
    quat2 = Quaternion(q2)
    assert quat1 * quat2 == Quaternion(expected)


@pytest.mark.parametrize("q, scalar, expected",
                         [[(-1, 3, 1, 2), 3, (-3, 9, 3, 6)],
                          [(-4, 0, 1, 5), 0, (0, 0, 0, 0)],
                          [(8, 13, 12, -16), -1, (-8, -13, -12, 16)],
                          [(-4, 10, 8, 15), 0.5, (-2, 5, 4, 7.5)]])
def test_multiply_scalar(q, scalar, expected):
    quat = Quaternion(q)
    assert quat * scalar == Quaternion(expected)


@pytest.mark.parametrize("q1, q2, expected",
                         [[(4, 6, -3, 3), (-8, 0, 2, 2), (0.278, -0.056, -0.833, -0.444)],
                          [(-8, 0, 2, 2), (4, 6, -3, 3), (-0.286, 0.057, 0.857, -0.457)]])
def test_divide_quaternion(q1, q2, expected):
    quat1 = Quaternion(q1)
    quat2 = Quaternion(q2)

    result = (quat1 / quat2).round(3)
    assert result == Quaternion(expected)


@pytest.mark.parametrize("q, scalar, expected",
                         [[(-1, 3, 1, 2), 3, (-0.3333, 1, 0.3333, 0.6667)],
                          [(8, 13, 12, -16), -1, (-8, -13, -12, 16)],
                          [(-4, 10, 8, 15), 2, (-2, 5, 4, 7.5)]])
def test_divide_scalar(q, scalar, expected):
    quat = Quaternion(q)
    assert (quat / scalar).round(4) == Quaternion(expected)


@pytest.mark.parametrize("q1, q2, expected", [[(1, 2, 3, 4), (5, 6, 7, 8), (6, 8, 10, 12)],
                                              [(1, 0, 1, 0), (0, 1, 0, 1), (1, 1, 1, 1)]])
def test_add_quaternion(q1, q2, expected):
    q1 = Quaternion(q1)
    q2 = Quaternion(q2)
    assert q1 + q2 == Quaternion(expected)


def test_add_scalar():
    q = Quaternion((1, 2, 3, 4))
    with pytest.raises(TypeError):
        q + 3


@pytest.mark.parametrize("q1, q2, expected", [[(1, 2, 3, 4), (5, 6, 7, 8), (-4, -4, -4, -4)],
                                              [(1, 0, 1, 0), (0, 1, 0, 1), (1, -1, 1, -1)]])
def test_subtract_quaternion(q1, q2, expected):
    q1 = Quaternion(q1)
    q2 = Quaternion(q2)
    assert q1 - q2 == Quaternion(expected)


def test_subtract_scalar():
    q = Quaternion((1, 2, 3, 4))
    with pytest.raises(TypeError):
        q - 3


def test_numpy_functions():
    quat = Quaternion((1, 2, 3, 4))
    array = np.asarray(quat)
    assert np.all(quat.scalar_last == array)
    assert quat.isscalar
    assert quat.shape == ()
    with pytest.raises(TypeError):
        len(quat)


def test_string_function():
    q = Quaternion((1, 2, 3, 4))
    assert str(q) == '(x, y, z, w) [1., 2., 3., 4.]'


#
# Test array version of quaternions
#

@pytest.mark.parametrize("scalar_last, x, y, z, w", [[[(1, 2, 3, 4), (5, 6, 7, 8)], (1, 5), (2, 6), (3, 7), (4, 8)]])
def test_component_properties_array(scalar_last, x, y, z, w):
    q = Quaternion(scalar_last)
    assert np.all(q.x == x)
    assert np.all(q.y == y)
    assert np.all(q.z == z)
    assert np.all(q.w == w)


@pytest.mark.parametrize("v1, v2, expected", [[[(1, 2, 3), (7, 8, 9)], [(4, 5, 6), (10, 11, 12)],
                                               [(-3, 6, -3, 64.83291), (-3, 6, -3, 532.10148)]]])
def test_from_vectors_array(v1, v2, expected):
    q = Quaternion.from_vectors(v1, v2)
    assert q.round(5) == Quaternion(expected)


def test_array_access():
    data = [(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12)]
    q = Quaternion(data)
    assert len(q) == 3
    assert q.shape == (3,)
    for i in range(3):
        assert q[i] == Quaternion(data[i])

    assert q[1:] == Quaternion(data[1:])
    assert q == Quaternion(data)
    assert q != Quaternion(data[1:])


def test_interpolation():
    quaternions = Quaternion([(0.4180201893760801, 0.8816012924578123, -0.016996745366792788, -0.21851634505829906),
                              (0.41793188408876947, 0.8815792563402776, -0.017441117860410714, -0.2187390282055626),
                              (0.4178444613661411, 0.881556498122404, -0.017885742276816587, -0.21896174767984167),
                              (0.41775611625340536, 0.8815340055270372, -0.01833111929389265, -0.21918392846593968),
                              (0.417667805843261, 0.8815112818574881, -0.018776785465131375, -0.21940577988658538)])

    obstimes = Time(['2022-07-27T23:59:00.740', '2022-07-27T23:59:01.740', '2022-07-27T23:59:02.740',
                     '2022-07-27T23:59:03.740', '2022-07-27T23:59:04.740'], format='isot', scale='utc')


def test_quaternion_attribute_quaternion():
    attr = QuaternionAttribute()

    q = Quaternion((0, 0, 0, 1))
    value, converted = attr.convert_input(q)
    assert value == q
    assert converted is False


def test_quaternion_attribute_array():
    attr = QuaternionAttribute()

    in_val = (0, 0, 0, 1)
    value, converted = attr.convert_input(in_val)
    assert isinstance(value, Quaternion)
    assert np.all(value.scalar_last == in_val)
    assert converted is True

    in_val = [(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12)]
    value, converted = attr.convert_input(in_val)
    assert isinstance(value, Quaternion)
    assert np.all(value.scalar_last == in_val)
    assert converted is True


def test_quaternion_attribute_invalid_array_size():
    attr = QuaternionAttribute()

    in_val = (1, 2, 3)
    with pytest.raises(TypeError):
        attr.convert_input(in_val)


def test_quaternion_attribute_invalid_value():
    attr = QuaternionAttribute()

    in_val = Time('2022-08-29T09:00:00', format='isot', scale='utc')
    with pytest.raises(TypeError):
        attr.convert_input(in_val)


def test_quaternion_attribute_none():
    attr = QuaternionAttribute()

    in_val = None
    value, converted = attr.convert_input(in_val)
    assert value is None
    assert converted is False
