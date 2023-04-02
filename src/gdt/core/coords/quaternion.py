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
from collections.abc import Sequence
from typing import Union, List
import numpy as np
from astropy.coordinates import Attribute
from astropy.time import Time
from scipy.spatial.transform import Rotation, Slerp
from astropy.utils.data_info import MixinInfo

from gdt.core.types import ArrayBase, assert_array_valid

VECTORS = Union[Sequence, np.ndarray]
QUATERNIONS = VECTORS
SCALARS = Union[int, float, Sequence, np.ndarray]


class Quaternion(ArrayBase):
    """ A class for containing a quaternion and performing quaternion operations.

    The quaternion is represented internally as a two-dimensional array of size 
    (n,4) where the elements follow numpy's row-major order (i.e. 
    array[row, colum]) and each quaternion is stored by row with the columns 
    representing i(x), j(y), k(z), and w in that order (scalar-last).
    
    Overloaded operations supported are:
    
    *  Addition
    *  Subtraction 
    *  Multiplication
    *  Division
    *  Negation
    *  Equivalence
    
    Parameters:
        quaternion (np.array): Either a 1D 4-element array containing a single
                               quaternion or a 2D array of shape (`n`, 4) 
                               containing `n` quaternions.
        scalar_first (bool, optional): Set to True if the input arrays are in
                                       scalar-first representation, False if in
                                       scalar-last representation.  Default is
                                       False.
    """
    # Attributes used by ArrayBase
    _rowsize = 4
    _prefix = '(x, y, z, w) '
    _dtype = np.float64
    _name = 'quaternion'

    info = MixinInfo()
    """(:class:`astropy.utils.data_info.MixinInfo`): Used by astropy Table to 
    store column information """


    EPSILON = 1e-12
    """(float): How close the floating point numbers need to be before they 
    are considered equal."""

    def __init__(self, quaternion: QUATERNIONS, scalar_first: bool = False):
        super().__init__(quaternion)
        if scalar_first:
            self._array = np.hstack((self._array[:, 1:4], self._array[:, 0].reshape(-1, 1)))

    @property
    def conjugate(self) -> 'Quaternion':
        """(:class:`Quaternion`): Conjugate of the quaternion(s)"""
        new = self._obj()
        new._array = np.hstack((-self._array[:, :3] + 0.0, self._array[:, 3].reshape(-1, 1)))
        return new

    @property
    def inverse(self) -> 'Quaternion':
        """(:class:`Quaternion`): The inverse of the quaternion(s)"""
        new = self.conjugate
        new._array = new._array / (self._array ** 2).sum(axis=1)[:, np.newaxis]
        return new

    @property
    def norm(self) -> Union[float, np.ndarray]:
        """(float or np.array): The norm of the quaternion(s)"""
        return self._scale(np.sqrt((self._array ** 2).sum(axis=1)))

    @property
    def rotation(self) -> Rotation:
        """(:class:`scipy.spatial.transform.Rotation`): The rotation matrix for 
        the quaternion(s)"""
        return Rotation(self.scalar_last)

    @property
    def scalar_first(self):
        """(np.array): Array of quaternion values with the scalar as the first 
        column"""
        val = np.hstack((self._array[:, 3].reshape(-1, 1), self._array[:, :3]))
        return self._scale(val)

    @property
    def scalar_last(self):
        """(np.array): An array of quaternion values with the scalar as 
        the last column"""
        return self._value

    @property
    def unit(self) -> 'Quaternion':
        """(:class:`Quaternion`): The quaternion(s) as unit quaternion(s)"""
        new = self._obj()
        new._array = self._array / np.asarray(self.norm).reshape(-1, 1)
        return new

    @property
    def w(self):
        """(np.array): The w value(s) in the quaternion"""
        return self._scale(self._array[:, 3])

    @property
    def xyz(self):
        """(np.array): The vector portion of the quaternion(s)"""
        return self._scale(self._array[:, :3])

    @property
    def x(self):
        """(np.array): The x value(s) in the quaternion"""
        return self._scale(self._array[:, 0])

    @property
    def y(self):
        """(np.array): The y value(s) in the quaternion"""
        return self._scale(self._array[:, 1])

    @property
    def z(self):
        """(np.array): The z value(s) in the quaternion"""
        return self._scale(self._array[:, 2])

    def dot(self, other):
        """Return the dot product of this quaternion and another.
        
        Args:
            other (:class:`Quaternion`): The other quaternion
        
        Returns:
            (:class:`Quaternion`)
        """
        return self._scale(np.einsum('ij,ij->i', self._array, other._array))

    def equal_rotation(self, other):
        """Determine if another quaternion represents the same rotation as this
        quaternion.
        
        Args:
            other (:class:`Quaternion`): The other quaternion
        
        Returns:
            (bool or np.array)
        """
    
        self_u = self.unit
        other_u = other.unit
        return np.all(np.abs(self_u.dot(other_u)) > 1.0 - self.EPSILON)

    def round(self, decimals: int) -> 'Quaternion':
        """Returns the quaternion(s) with the components rounded to the given decimal places

        Args:
            decimals: Number of decimal places to round the components of the quaternion.

       Returns:
            (:class:`Quaternion`)
         """
        new = self._obj()
        new._array = np.round(self._array, decimals=decimals)
        return new

    @classmethod
    def from_rotation(cls, rot: Rotation) -> 'Quaternion':
        """Create a quaternion from SciPy Rotation objects.

         Parameters:
            rot (scipy.spatial.transform.Rotation): The SciPy Rotation object

        Returns:
            (:class:`Quaternion`)
        """
        return cls(rot.as_quat(), scalar_first=False)

    @classmethod
    def from_vectors(cls, vec1: Union[np.ndarray, List], vec2: Union[np.ndarray, List]):
        """Create a quaternion from the rotation between two vectors.

        Args:
            vec1 (np.array): A cartesian vector
            vec2 (np.array): Another cartesian vector
        
        Returns:
            (:class:`Quaternion`)
        """
        # verify dimensions and type
        vec1 = np.asarray(vec1)
        vec2 = np.asarray(vec2)
        assert_array_valid(vec1, rowsize=3, name='vec1')
        assert_array_valid(vec2, rowsize=3, name='vec2')
        if vec1.shape == (3,):
            vals = np.hstack((np.cross(vec1, vec2),
                              np.linalg.norm(vec1) * np.linalg.norm(vec2) + np.dot(vec1, vec2)))
        else:
            vals = np.hstack((np.cross(vec1, vec2),
                              np.asarray([np.linalg.norm(vec1, axis=1) * np.linalg.norm(vec2, axis=1) +
                                          np.einsum('ij,ij->i', vec1, vec2)]).T))
        obj = cls(vals)
        return obj

    @classmethod
    def from_xyz_w(cls, xyz: VECTORS, w: SCALARS) -> 'Quaternion':
        """Create a quaternion from a vector and scalar.

        Args:
            xyz: Represents the values for the x, y and z axis
            w: Represents the rotation for the resultant (w)

        Returns:
           (:class:`Quaternion`)
        """
        xyz = np.asarray(xyz)
        w = np.asarray(w)
        assert_array_valid(xyz, rowsize=3, name='xyz')
        assert_array_valid(w, rowsize=1, name='w')
        if xyz.shape == (3,) and w.shape == ():
            value = np.concatenate((xyz, w.reshape(1,)))
        else:
            value = np.hstack((xyz.reshape(-1, 3), w.reshape(-1, 1)))
        return cls(value)

    def _setup_interpolation(self, unix_tai: Union[float, np.ndarray]) -> Slerp:
        """Called to initialize the interpolation function"""
        return Slerp(unix_tai, Rotation.from_quat(self))

    def __add__(self, other):
        """Perform addition.
        Only the addition of another Quaternion is allowed.

        Parameters:
            other: The addend
        """
        if not isinstance(other, self.__class__):
            raise TypeError('Only addition of two quaternions are supported.')
        new = self._obj()
        new._array = self._array + other._array
        return new

    def __eq__(self, other):
        """Determine if the two quaternions equate to the same rotation."""
        if self._array.shape != other._array.shape:
            return False
        return np.allclose(self._array, other._array, atol=self.EPSILON)

    def __mul__(self, other):
        """Perform multiplication.
        If it's Quaternion * Quaternion, then Hamilton Product is used to speed up the operation.
        Otherwise, scalar multiplication is performed following Numpy broadcast rules.

        Parameters:
            other: The multiplier
        """
        new = self._obj()
        if isinstance(other, self.__class__):
            s = self._array
            o = other._array
            new._array = np.hstack((
                # Using Hamilton Product to reduce the time needed by half
                # The resulting stack will be in (x,y,z,w) order with the components having the following indexes:
                # x = [:,0], y = [:,1], z = [:,2], w = [:,3]
                (s[:, 3] * o[:, 0] + s[:, 0] * o[:, 3] + s[:, 1] * o[:, 2] - s[:, 2] * o[:, 1])[:, np.newaxis],
                (s[:, 3] * o[:, 1] - s[:, 0] * o[:, 2] + s[:, 1] * o[:, 3] + s[:, 2] * o[:, 0])[:, np.newaxis],
                (s[:, 3] * o[:, 2] + s[:, 0] * o[:, 1] - s[:, 1] * o[:, 0] + s[:, 2] * o[:, 3])[:, np.newaxis],
                (s[:, 3] * o[:, 3] - s[:, 0] * o[:, 0] - s[:, 1] * o[:, 1] - s[:, 2] * o[:, 2])[:, np.newaxis]
            ))
        else:
            new._array = self._array * other
        return new

    def __neg__(self):
        """Returns the negative of the quaternion(s)"""
        new = self._obj()
        new._array = -self._array
        return new

    def __sub__(self, other):
        """Perform subtraction.
        Only the subtraction of another Quaternion is allowed.

        Parameters:
            other: The subtrahend
        """
        if not isinstance(other, self.__class__):
            raise TypeError('Only addition of two quaternions are supported.')
        new = self._obj()
        new._array = self._array - other._array
        return new

    def __truediv__(self, other) -> 'Quaternion':
        """Perform division.
        If it's Quaternion / Quaternion, then the operation is performed by multiplying with the inverse of other(s).
        Otherwise, scalar division is performed following Numpy broadcast rules.

        Parameters:
            other: The divisor
        """
        if isinstance(other, self.__class__):
            inv = other.inverse
            new = self.__mul__(inv)
        else:
            new = self._obj()
            new._array = self._array / other
        return new


class QuaternionAttribute(Attribute):
    """
    A quaternion attribute class to use with
    astropy.coordinates.BaseCoordinateFrame
    """

    def convert_input(self, value):
        """Function called by Astropy Representation/Frame system.
        This function verifies that the value is a quaternion.
        """
        if value is None:
            converted = False
        elif isinstance(value, (tuple, list, np.ndarray)):
            try:
                value = Quaternion(value)
                converted = True
            except (ValueError, TypeError):
                raise TypeError('Value must be a Quaternion object or an array that can be converted to a quaternion.')
        elif isinstance(value, Quaternion):
            converted = False
        else:
            raise TypeError('Value must be a Quaternion object or an array that can be converted to a quaternion.')

        return value, converted
