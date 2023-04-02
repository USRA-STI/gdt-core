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
import os
import numpy as np

__all__ = ['SouthAtlanticAnomaly']

class SouthAtlanticAnomaly():
    """A base class for a South Atlantic Anomaly (SAA) boundary in Earth 
    latitude and longitude. This class should be further sub-classed with the 
    following the class variables:
    
        *  _latitude - A list or array of latitude points
        *  _longitude - A list or array of longitude points    
    """
    def __init__(self):   
        # _latitude and _longitude must be set
        if not hasattr(self, '_latitude'):
            raise AttributeError("SouthAtlanticAnomaly must have class " \
                                 "attribute '_latitude' defined")
        if not hasattr(self, '_longitude'):
            raise AttributeError("SouthAtlanticAnomaly must have class " \
                                 "attribute '_longitude' defined")
    
        self._latitude = np.asarray(self._latitude)
        self._longitude = np.asarray(self._longitude)
        if self._latitude.size != self._longitude.size:
            raise ValueError('_longitude and _latitude must be the same size')
    
    @property
    def latitude(self):
        """(np.array): The latitude points of the boundary"""
        return self._latitude
    
    @property
    def longitude(self):
        """(np.array): The longitude points of the boundary"""
        return self._longitude
    
    @property
    def num_points(self):
        """(int): Number of vertices in polygon"""
        return self.latitude.size
    
    def is_closed(self):
        """Test if the boundary represents a closed polygon.  A closed polygon
        is one where the first and last points are equal.
        
        Returns:
            (bool)
        """
        if (self.latitude[0] == self.latitude[-1]) & \
           (self.longitude[0] == self.longitude[-1]):
            return True
        else:
            return False
    
    def is_convex(self):
        """Test if the boundary represents a convex polygon. The test is 
        performed by measuring every interior angle of the polygon and checking
        if all angles are < 180 deg.  The definition of a convex polygon is one
        where all interior angles are < 180 deg.
        
        Returns:
            (bool)
        """
        # calculate the vectors
        dlat = self.latitude[1:] - self.latitude[:-1]
        dlon = self.longitude[1:] - self.longitude[:-1]
        dz = np.zeros_like(dlon)
        vec = np.vstack((dlon, dlat, dz))
        
        norms = np.linalg.norm(vec, axis=0)
        
        winding = 0.0
        angles = np.zeros(dlat.size-1)
        for i in range(dlat.size-1):
            sign = np.sign(np.cross(vec[:,i], vec[:,i+1])[-1])
            angle = sign * np.arccos(np.dot(vec[:,i], vec[:,i+1]) / \
                                     (norms[i] * norms[i+1]))
            angles[i] = angle
            winding += angle
        
        if winding < 0.0:
            angles += np.pi
        else:
            angles = np.pi - angles
        
        if (angles >= np.pi).sum() > 0:
            return False
        else:
            return True
