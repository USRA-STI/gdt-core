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
from enum import Enum
import astropy.units as u

__all__ = ['Detectors']

from astropy.coordinates import SkyCoord

# unfortunately Sphinx has a major bug that prevents the autodoc of Enums, 
# so we have to define all of this in the docstring...

class Detectors(Enum):
    """Detector base class containing the name/number mapping and the 
    pointing direction in spacecraft coordinates.
    
    When inheriting this base class, each detector will be a 4-tuple defined
    as a class variable, with the class variable name represented as the 
    short name of the detector.  The 4-tuple values are the full detector name,
    the corresponding detector number, and the detector azimuth and zenith.
    
    Example:
    
        >>> class MyDetectors(Detectors):
        >>>     mydet0 = ('MyDetector0', 0, 10.0, 50.0)
        >>>     mydet1 = ('MyDetector1', 1, 0.0, 30.0)
    
    Parameters:
        full_name (str): The full detector name
        number (int): The detector number
        azimuth (float): The azimuth of the detector normal
        zenith (float): The zenith of the detector normal
    
    .. rubric:: Attributes Summary
    .. autosummary::

      azimuth
      elevation
      full_name
      number
      zenith

    .. rubric:: Methods Summary

    .. autosummary::

      from_full_name
      from_num
      from_str
      pointing
      skycoord
  
    .. rubric:: Attributes Documentation

    .. autoattribute:: azimuth
    .. autoattribute:: elevation
    .. autoattribute:: full_name
    .. autoattribute:: number
    .. autoattribute:: zenith

    .. rubric:: Methods Documentation

    .. automethod:: from_full_name
    .. automethod:: from_num
    .. automethod:: from_str
    .. automethod:: pointing
    .. automethod:: skycoord
    """
    def __init__(self, full_name, number, azimuth, zenith):
        self._full_name = full_name
        self._number = number
        self._azimuth = azimuth
        self._zenith = zenith

    @property
    def azimuth(self):
        """(float): The azimuth of the detector normal"""
        return self._azimuth

    @property
    def elevation(self):
        """(float): The elevation of the detector normal"""
        return 90 * u.deg - self._zenith

    @property
    def full_name(self):
        """(str): The full detector name"""
        return self._full_name

    @property
    def number(self):
        """(int): The detector number"""
        return self._number

    @property
    def zenith(self):
        """(float): The zenith of the detector normal"""
        return self._zenith

    def pointing(self):
        """The detector pointing in azimuth and zenith
        
        Returns:
            (float, float)
        """
        return self.azimuth, self.zenith

    def skycoord(self, frame):
        """Return an Astropy SkyCoord of the detector pointing in the
        SpacecraftFrame.
        
        Args:
            frame (:class:`~gdt.core.coords.SpacecraftFrame`): 
                The spacecraft frame
        
        Returns:
            (astropy.coordinates.SkyCoord)
        """
        return SkyCoord(self.azimuth, self.elevation, frame=frame)

    @classmethod
    def from_full_name(cls, full_name):
        """Create a Detector from the full detector name
        
        Args:
            full_name (str): The full name of the detector
        
        Returns:
            (:class:`Detector`)
        """
        for d in cls:
            if d.full_name == full_name:
                return d
        raise ValueError('{} not a valid detector'.format(full_name))

    @classmethod
    def from_num(cls, num):
        """Create a Detector from an index number
        
        Args:
            num (int): The index number
        
        Returns:
            (:class:`Detector`)
        """
        for d in cls:
            if d.number == num:
                return d
        raise ValueError('{} not a valid detector number'.format(num))

    @classmethod
    def from_str(cls, value):
        """Create a Detector from a short string name
        
        Args:
            value (str): The detector name
        
        Returns:
            (:class:`Detector`)
        """
        if value in cls.__members__:
            return cls[value]
        raise ValueError('{} not a valid detector string'.format(value))
    
    def __repr__(self):
        return "<{0}: {1}>".format(self.__class__.__name__, self.name)
