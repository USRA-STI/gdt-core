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
from typing import Optional

import numpy as np
from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME
import astropy.coordinates as a_coords
import astropy.coordinates.representation as r
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.timeseries import TimeSeries
from scipy.interpolate import interp1d
from scipy.spatial.transform import Slerp, Rotation
from astropy.constants import R_earth

from ..quaternion import QuaternionAttribute, Quaternion
from gdt.core.types import as_times
from gdt.core.detector import Detectors


class DetectorsAttribute(a_coords.Attribute):
    """
    A detectors attribute class to use with
    astropy.coordinates.BaseCoordinateFrame
    """

    def convert_input(self, value):
        """Function called by Astropy Representation/Frame system.
        This function verifies that the value is a Detectors object.
        """
        if value is None:
            converted = False
        elif issubclass(value, Detectors):
            converted = False
        else:
            raise TypeError('Value must be a Detectors object.')

        return value, converted


class SpacecraftFrame(a_coords.BaseCoordinateFrame):
    """
    Spacecraft frame. The frame is defined by spacecraft's quaternion that 
    represents a rotation from the spacecraft frame to the ICRS frame as well
    as the spacecraft's position in Earth-centered Inertial Coordinates.
    
    Parameters:
        obstime (astropy.time.Time): The time(s) at which the spacecraft frame
                                     is represented.
        obsgeoloc (astropy.coordinates.representation.CartesianRepresentation, optional):
            The position of the spacecraft in orbit in Earth-centered Inertial
            coordinates.
        obsgeovel (astropy.coordinates.representation.CartesianRepresentation, optional):
            The orbital velocity of the spacecraft in Earth-centered Inertial
            coordinates.
        quaternion (:class:`Quaternion`, optional): 
            The attitude quaternions of the spacecraft representing the 
            spacecraft orientation. 
        detectors (:class:`~gdt.core.detector.Detectors`, optional): 
            Detector definitions associated with this frame.
    """
    default_representation = r.SphericalRepresentation
    """The default representation of the spacecraft frame"""

    frame_specific_representation_info = {
        "spherical": [a_coords.RepresentationMapping("lon", "az", defaultunit=u.degree),
                      a_coords.RepresentationMapping("lat", "el", defaultunit=u.degree)]
    }

    obstime = a_coords.TimeAttribute(default=DEFAULT_OBSTIME)
    """Time or times of associated with the spacecraft frame"""
    
    obsgeoloc = a_coords.CartesianRepresentationAttribute(default=[0, 0, 0],
                                                          unit=u.m)
    """Spacecraft position in Earth-centered Inertial coordinates at time(s)
    ``obstime``"""
                
    obsgeovel = a_coords.CartesianRepresentationAttribute(default=[0, 0, 0],
                                                          unit=u.m / u.s)
    """Spacecraft velocity at time(s) ``obstime``"""
    
    quaternion = QuaternionAttribute(default=None)
    """Spacecraft attitude quaternions at time(s) ``obstime``"""
    
    detectors = DetectorsAttribute(default=None)
    """The spacecraft detectors definition"""

    def __init__(self, *args, copy=True, representation_type=None, 
                 differential_type=None, **kwargs):
        super().__init__(*args, copy=copy, 
                         representation_type=representation_type, 
                         differential_type=differential_type,
                         **kwargs)

        # some background interpolation
        self._interp_geoloc: Optional[interp1d] = None
        self._interp_geovel: Optional[interp1d] = None
        self._interp_quat: Optional[Slerp] = None
        self._interp: bool = False

    @property
    def earth_angular_radius(self):
        """The apparent angular radius, in degrees, of the Earth from the
        spacecraft position

        Returns:
            (astropy.units.Quantity)
        """
        return np.rad2deg(np.arcsin(R_earth / (R_earth + \
                                               self.earth_location.height)))

    @property
    def earth_location(self) -> a_coords.EarthLocation:
        """(astropy.coordinates.EarthLocation): The spacecraft position as an
        Earth coordinate."""
        return a_coords.EarthLocation(*self.sc_itrs.cartesian.xyz)

    @property
    def geocenter(self) -> a_coords.GCRS:
        """(astropy.coordinates.GCRS): The geocenter location in the GCRS 
        frame."""
        return a_coords.GCRS(-self.obsgeoloc, obstime=self.obstime)

    @property
    def sc_gcrs(self) -> a_coords.GCRS:
        """(astropy.coordinates.GCRS): The spacecraft position in the GCRS 
        frame"""
        return a_coords.GCRS(self.obsgeoloc, obstime=self.obstime)

    @property
    def sc_itrs(self) -> a_coords.ITRS:
        """(astropy.coordinates.ITRS): The spacecraft position in the ITRS 
        frame"""
        return self.sc_gcrs.transform_to(a_coords.ITRS(obstime=self.obstime))

    @property
    def sun_visible(self):
        """(bool): True if the sun is visible to the spacecraft."""
        return self.location_visible(a_coords.get_sun(self.obstime))

    def at(self, obstime: Time):
        """Retrieve the interpolated spacecraft positions and quaternions for 
        the specified time(s).
        
        Args:
            obstime (astropy.time.Time): The times for which the frames are 
                                         requested.
        
        Returns:
            (:class:`SpacecraftFrame`)
        """
        if not self._interp:
            self.init_interpolation()

        if self._interp_geoloc is None:
            raise ValueError('Interpolation can not be performed on obsgeoloc.')

        t = obstime.unix_tai

        if self._interp_geovel is None:
            geovel = None
        else:
            geovel = r.CartesianRepresentation(self._interp_geovel(t), unit=self.obsgeovel.x.unit)

        if self._interp_quat is None:
            quat = None
        else:
            quat = Quaternion.from_rotation(self._interp_quat(t))

        obj = self.__class__(
            obstime=obstime,
            obsgeoloc=r.CartesianRepresentation(self._interp_geoloc(t), unit=self.obsgeoloc.x.unit),
            obsgeovel=geovel,
            quaternion=quat,
            detectors=self.detectors
        )
        return obj

    def detector_angle(self, det, coord):
        """Calculate the detector angle for a SkyCoord.
        
        Note:
            This will only work if the SpacecraftFrame was initialized with a
            valid :class:`~gdt.core.detector.Detectors` object containing the 
            detector pointing definitions.
        
        Args:
            det (string or :class:`~gdt.core.detector.Detectors`):
                The detector for which to calculate the angle.  If this is a
                string value, it can either be the ``name`` or ``full_name`` of
                the detector.  A single Detector object can also be used.
            coord (astropy.coordinates.SkyCoord): 
                The coordinate(s) against which the angle will be calculated.
        
        Returns:
            (astropy.coordinates.Angle)
        """
        if self.detectors is None:
            raise AttributeError('The detectors attribute was not set on ' \
                                 'initialization')
        if isinstance(det, self.detectors):
            det_obj = det
        else:
            # try the simple name first. if not found, try the full name
            try:
                det_obj = self.detectors.from_str(det)
            except:
                det_obj = self.detectors.from_full_name(det)
        
        det_pointing = (det_obj.azimuth, det_obj.elevation)
        if self.ndim == 0:
            det_coord = a_coords.SkyCoord(*det_pointing, frame=self)
        else:
            det_coord = a_coords.SkyCoord([det_pointing] * self.shape[0], 
                                          frame=self)
        return det_coord.separation(coord)
    
    def estimate_velocity(self, time: Time, eps=0.1):
        """Estimate the orbital velocity at a given time. The result is returned
        as an astropy Quantity in units of `x`/s, where `x` is the distance unit
        of the Earth-inertial coordinates used to initialize the object.

        Args:
            time (float or astropy.time.Time): The time(s) at which to
                                               estimate the orbital velocity.
            eps (float): The epsilon change in time, in seconds, over which the
                         time derivative is calculated. Default is 0.1.

        Returns:
            (:class:`astropy.units.Quantity`)
        """
        if not isinstance(time, (Time, TimeSeries)):
            raise TypeError('time must be a astropy.Time object')
        input_scalar = time.isscalar

        time = as_times(time)
        tdelta = TimeDelta(eps, scale='tai', format='sec') / 2.0
        time1 = time - tdelta
        time2 = time + tdelta

        eic1 = self.at(time1).obsgeoloc
        eic2 = self.at(time2).obsgeoloc

        vel = (eic2 - eic1) / (eps * u.s)
        if input_scalar and len(vel) == 1:
            vel = vel[0]
        return vel

    def init_interpolation(self):
        """Initialize the interpolation functions"""
        t = self.obstime.unix_tai
        if self.obstime is not None and not self.obstime.isscalar and len(self.obstime) > 1:
            if self.quaternion is not None:
                self._interp_quat = Slerp(t, Rotation.from_quat(self.quaternion))
            if self.obsgeoloc is not None and self.obsgeoloc.shape != ():
                self._interp_geoloc = interp1d(t, self.obsgeoloc.xyz.value)
            if self.obsgeovel is not None and self.obsgeovel.shape != ():
                self._interp_geovel = interp1d(t, self.obsgeovel.xyz.value)

            self._interp = True

    def location_visible(self, pos: a_coords.SkyCoord):
        """Determine if a sky location is visible or occulted by the Earth.

        Args:
            pos (astropy.coordinates.SkyCoord): Sky location(s)

        Returns:
            np.array(dtype=bool): A boolean array where True indicates the 
            position is visible.
        """
        return self.geocenter.separation(pos) > self.earth_angular_radius

    def __repr__(self):
        try:
            sz = self.shape[0]
        except:
            sz = 1
        s = '<{0}: {1} frames;\n'.format(self.__class__.__name__, sz)
        if sz == 1:
            s += ' obstime=[{}]\n'.format(self.obstime.__str__())
            try:
                s += ' obsgeoloc=[{}]\n'.format(self.obsgeoloc.__str__())
            except:
                pass
            try:
                s += ' obsgeovel=[{}]\n'.format(self.obsgeovel.__str__())
            except:
                pass
            try:
                s += ' quaternion=[{}]>'.format(self.quaternion.__str__())
            except:
                s += '>'
        else:
            s += ' obstime=[{}, ...]\n'.format(self.obstime[0].__str__())
            try:
                s += ' obsgeoloc=[{}, ...]\n'.format(self.obsgeoloc[0].__str__())
            except:
                pass
            try:
                s += ' obsgeovel=[{}, ...]\n'.format(self.obsgeovel[0].__str__())
            except:
                pass
            try:
                s += ' quaternion=[{}, ...]>'.format(self.quaternion[0].__str__())
            except:
                s += '>'
        return s
        

@a_coords.frame_transform_graph.transform(a_coords.FunctionTransform, SpacecraftFrame, a_coords.ICRS)
def spacecraft_to_icrs(sc_frame, icrs_frame):
    """Convert from the spacecraft frame to the ICRS frame.

    Args:
        sc_frame (:class:`SpacecraftFrame`): The spacecraft frame
        icrs_frame (astropy.coordinates.ICRS)

    Returns:
        (astropy.coordinates.ICRS)
    """
    xyz = sc_frame.cartesian.xyz.value
    rot = Rotation.from_quat(sc_frame.quaternion)
    xyz_prime = rot.apply(xyz.T)
    if xyz_prime.ndim == 1:
        xyz_prime = xyz_prime.reshape(1, -1)
    ra = np.arctan2(xyz_prime[:, 1], xyz_prime[:, 0])
    mask = (ra < 0.0)
    ra[mask] += 2.0 * np.pi
    dec = np.pi / 2.0 - np.arccos(xyz_prime[:, 2])

    icrs = a_coords.ICRS(ra=ra * u.radian, dec=dec * u.radian)
    return icrs.transform_to(icrs_frame)


@a_coords.frame_transform_graph.transform(a_coords.FunctionTransform, a_coords.ICRS, SpacecraftFrame)
def icrs_to_spacecraft(icrs_frame, sc_frame):
    """Convert from the ICRS frame to the spacecraft frame.

    Args:
        icrs_frame (astropy.coordinates.ICRS)
        sc_frame (:class:`SpacecraftFrame`): The spacecraft frame

    Returns:
        (:class:`SpacecraftFrame`)
    """
    xyz = icrs_frame.cartesian.xyz.value
    rot = Rotation.from_quat(sc_frame.quaternion)
    xyz_prime = rot.inv().apply(xyz.T)
    if xyz_prime.ndim == 1:
        xyz_prime = xyz_prime.reshape(1, -1)
    az = np.arctan2(xyz_prime[:, 1], xyz_prime[:, 0])
    mask = (az < 0.0)
    az[mask] += 2.0 * np.pi
    el = np.pi / 2.0 - np.arccos(xyz_prime[:, 2])
    return type(sc_frame)(az=az * u.radian, el=el * u.radian,
                          quaternion=sc_frame.quaternion)
