# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract Nos.: CA 80MSFC17M0022 / 80NSSC24M0035
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
import numpy as np

__all__ = ['Parameter']


class Parameter():
    """A fit parameter class
    
    Parameters:
        value (float): The central fit value
        uncert (float or 2-tuple): The 1-sigma uncertainty. If a 2-tuple, then 
                                   is of the form (low, high)
        name (str, optional): The name of the parameter
        units (str, optional): The units of the parameter
        support (2-tuple, optional): The valid support of the parameter
    """
    def __init__(self, value, uncert, name='', units=None,
                 support=(-np.inf, np.inf)):

        self._value = float(value)
        if isinstance(uncert, (tuple, list)):
            if len(uncert) == 2:
                pass
            elif len(uncert) == 1:
                uncert = (uncert[0], uncert[0])
            else:
                raise ValueError('uncertainty must be a 1- or 2-tuple')
        elif isinstance(uncert, float):
            uncert = (uncert, uncert)
        else:
            raise TypeError('uncertainty must be a float or 1- or 2-tuple')
        self._uncert = uncert

        self._units = units
        self._name = name
        self._support = support

    @property
    def name(self):
        """(str): The name of the parameter"""
        return self._name

    @property
    def support(self):
        """(2-tuple): The valid support of the parameter"""
        return self._support

    @property
    def uncertainty(self):
        """(float, float): The 1-sigma uncertainty"""
        return self._uncert

    @property
    def units(self):
        """(str): The units of the parameter"""
        return self._units

    @property
    def value(self):
        """(float): The central fit value"""
        return self._value

    def one_sigma_range(self):
        """Return the 1 sigma range of the parameter fit
        
        Returns:
            (tuple): 2-tuple (low, high)
        """
        return (self.value - self.uncertainty[0],
                self.value + self.uncertainty[1])

    def to_fits_value(self):
        """Return as a tuple to be used for a FITS file

        Returns:
            (tuple): 3-value tuple (value, +uncertainty, -uncertainty)   
        """
        return (self.value, *self.uncertainty[::-1])

    def valid_value(self):
        """Check if the parameter value is within the allowed parameter range
        
        Returns:
            (bool)
        """
        if (self.value >= self.support[0]) and \
                (self.value <= self.support[1]):
            return True
        else:
            return False

    def _str_format(self):
        if (self.value > 0.005) and (self.uncertainty[0] > 0.005):
            value = '{0:.2f}'.format(self.value)
            uncertainty = tuple(
                ['{0:.2f}'.format(u) for u in self.uncertainty])
        else:
            value = '{0:.2e}'.format(self.value)
            val_coeff, val_exp = value.split('e')
            val_exp = int(val_exp)
            uncertainty = ['{0:.2e}'.format(u) for u in self.uncertainty]
            uncert_coeff = []
            uncert_exp = []
            for uncert in uncertainty:
                uncert_coeff.append(uncert.split('e')[0])
                uncert_exp.append(int(uncert.split('e')[1]))
        return (value, uncertainty)

    def __repr__(self):
        return '<{0}: {1}>'.format(self.__class__.__name__, self.name)

    def __str__(self):
        value, uncertainty = self._str_format()

        if uncertainty[0] == uncertainty[1]:
            s = '+/- {0}'.format(uncertainty[0])
        else:
            s = '+{1}/-{0}'.format(uncertainty[0], uncertainty[1])
        if self.units is None:
            return '{0}: {1} {2}'.format(self.name, value, s)
        else:
            return '{0}: {1} {2} {3}'.format(self.name, value, s, self.units)
