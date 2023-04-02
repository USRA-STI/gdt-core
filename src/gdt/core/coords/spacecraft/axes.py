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
import numpy as np
from astropy.coordinates import SkyCoord, Attribute


class SpacecraftAxes:
    """
    A class containing the axes pointings of a spacecraft in ICRS coordinates

    Parameters:
        x_pointing (SkyCoord, optional): The x-axis pointing
        y_pointing (SkyCoord, optional): The y-axis pointing
        z_pointing (SkyCoord, optional): The z-axis pointing
    """

    def __init__(self, x_pointing=None, y_pointing=None, z_pointing=None):

        if isinstance(x_pointing, SkyCoord) or x_pointing is None:
            self._x = x_pointing
        else:
            raise TypeError('x_pointing must be a SkyCoord or None')

        if isinstance(y_pointing, SkyCoord) or y_pointing is None:
            self._y = y_pointing
        else:
            raise TypeError('y_pointing must be a SkyCoord or None')

        if isinstance(z_pointing, SkyCoord) or z_pointing is None:
            self._z = z_pointing
        else:
            raise TypeError('z_pointing must be a SkyCoord or None')

    @property
    def has_x(self):
        """(bool): True if the x_pointing is set"""
        return True if self._x is not None else False

    @property
    def has_y(self):
        """(bool): True if the y_pointing is set"""
        return True if self._y is not None else False

    @property
    def has_z(self):
        """(bool): True if the z_pointing is set"""
        return True if self._z is not None else False

    @property
    def x_vector(self):
        """(np.array): The x-axis reference unit vector"""
        return np.array([1.0, 0.0, 0.0])

    @property
    def y_vector(self):
        """(np.array): The x-axis reference unit vector"""
        return np.array([0.0, 1.0, 0.0])

    @property
    def z_vector(self):
        """(np.array): The x-axis reference unit vector"""
        return np.array([0.0, 0.0, 1.0])

    def pointing_vector(self, axis):
        """The pointing vector of the requested axis

        Args:
            axis (str): Either 'x', 'y', or 'z'

        Returns:
            (np.array)
        """
        if axis == 'x':
            if self._x is not None:
                return self._x.cartesian.xyz.value
        elif axis == 'y':
            if self._y is not None:
                return self._y.cartesian.xyz.value
        elif axis == 'z':
            if self._z is not None:
                return self._z.cartesian.xyz.value
        else:
            raise ValueError('axis must be either x, y, or z')

    def __repr__(self):  # pragma: no cover - Repr is only used for debugging.
        x_str = ''
        y_str = ''
        z_str = ''
        if self.has_x:
            x_str = 'X axis: RA {0}, Dec {1}\n'.format(self._x.ra.value,
                                                       self._x.dec.value)
        if self.has_y:
            y_str = 'Y axis: RA {0}, Dec {1}\n'.format(self._y.ra.value,
                                                       self._y.dec.value)
        if self.has_z:
            z_str = 'Z axis: RA {0}, Dec {1}'.format(self._z.ra.value,

                                                     self._z.dec.value)
        s = '<SpacecraftAxes>\n' + x_str + y_str + z_str

        return s


class SpacecraftAxesAttribute(Attribute):
    """
    A spacecraft axes attribute class to use with
    astropy.coordinates.BaseCoordinateFrame
    """

    def convert_input(self, value):
        """Function called by Astropy Representation/Frame system.
        This function verifies that the value is a SpacecraftAxes.
        """
        if value is None:
            return None, False

        if not isinstance(value, SpacecraftAxes):
            raise TypeError('value must be a SpacecraftAxes object')

        converted = False
        return value, converted
