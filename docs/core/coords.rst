.. _core-coords:
.. |Quaternion| replace:: :class:`~gdt.core.coords.Quaternion`
.. |SpacecraftAxes| replace:: :class:`~gdt.core.coords.SpacecraftAxes`
.. |SpacecraftFrame| replace:: :class:`~gdt.core.coords.SpacecraftFrame`
.. |SpacecraftPosition| replace:: :class:`~gdt.core.coords.SpacecraftPosition`
.. |FitsFileContextManager| replace:: :class:`~gdt.core.file.FitsFileContextManager`
.. |SpacecraftFrameModelMixin| replace:: :class:`~gdt.core.coords.spacecraft.SpacecraftFrameModelMixin`
.. |SpacecraftStatesModelMixin| replace:: :class:`~gdt.core.coords.spacecraft.SpacecraftStatesModelMixin`
.. |BaseCoordinateFrame| replace:: BaseCoordinateFrame
.. _BaseCoordinateFrame: https://docs.astropy.org/en/stable/api/astropy.coordinates.BaseCoordinateFrame.html

***********************************************************************
Spacecraft Attitude, Position, and Coordinates (:mod:`gdt.core.coords`)
***********************************************************************

Introduction
============
Auxiliary mission data such as spacecraft/instrument attitude and position at
a given time are important pieces of information for a wide variety of reasons,
including determination if a source is visible, and which detector(s) should be
used to analyze a source.  This section introduces the basic classes that
spacecraft/instruments can use to interface with spacecraft attitude and 
position history data. 

Quaternions
===========
Often spacecraft provide attitude information in the form of a *quaternion*, 
a 4-element representation of a 3-dimensional rotation containing a vector part
(typically notated as `x`, `y`, `z`) and a scalar part (`w`). Quaternions are an 
information-efficient way to describe rotations, or series of rotations, that 
are not susceptible to gimbal lock. The GDT provides a |Quaternion| class to 
perform quaternion operations and use quaternions to define a spacecraft 
coordinate frame that is discussed below.

Examples
--------

There are a few different ways to create a |Quaternion| object.  One way is to
initialize with a 4-element vector representing the quaternion:

    >>> from gdt.core.coords import Quaternion
    >>> quat = Quaternion([0.0, 0.0, 1.0, 1.0], scalar_first=False)
    >>> quat
    <Quaternion (x, y, z, w)  [0., 0., 1., 1.] >

Notice that we initialize the object with an optional keyword 
``scalar_first=False.`` There are two conventions for representing quaternions: 
with the scalar part listed first, or with the scalar part listed last. 
To keep track of the where the scalar part is listed, we set the 
``scalar_first`` flag to True or False. 

You can also create a Quaternion object by explicitly specifying the vector and
scalar parts:

    >>> quat = Quaternion.from_xyz_w(xyz=[0.0, 0.0, 1.0], w=[1.0])
    >>> quat
    <Quaternion (x, y, z, w)  [[0., 0., 1., 1.]] >

A Quaternion object can also be created from the rotation between two Cartesian
vectors:

    >>> Quaternion.from_vectors([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    <Quaternion (x, y, z, w)  [0., 0., 1., 1.] >
   
There are a variety of properties and functions we can apply:

      >>> # scalar and vector parts
      >>> quat.w
      array([1.])
      >>> quat.xyz
      array([[0., 0., 1.]])
  
      >>> # the normalization of the quaternion
      >>> quat.norm
      array([1.41421356])
  
      >>> # the conjugate
      >>> quat.conjugate
      <Quaternion (x, y, z, w)  [[ 0.,  0., -1.,  1.]] >
  
      >>> # the inverse
      >>> quat.inverse()
      <Quaternion (x, y, z, w)  [[ 0. ,  0. , -0.5,  0.5]] >
  
      >>> # convert to a scalar-first or scalar-last representation
      >>> quat.scalar_first
      array([[1., 0., 0., 1.]])
      >>> quat.scalar_last
      array([[0., 0., 1., 1.]])  
  
      >>> # normalize and return unit quaternion
      >>> quat.unit
      <Quaternion (x, y, z, w)  [[0., 0., 0.70710678, 0.70710678]] >
  

Quaternion multiplication is a useful operation for describing a series of 
rotations.  You can simply apply quaternion rotations with the overloaded 
``*`` operator:

    >>> quat2 = Quaternion([1.0, 0.0, 1.0, 0.0], scalar_first=True)
    >>> quat * quat2
    <Quaternion (x, y, z, w)  [[-1.,  1.,  1.,  1.]] >
    >>> quat2 * quat
    <Quaternion (x, y, z, w)  [1., 1., 1., 1.] >

Notice that quaternion multiplication is not commutative, so be careful to keep
track of the order of rotations. Other overloaded operations include adding, subtracting, and dividing 
quaternions.  Quaternions can be determined to be either directly equivalent:

    >>> quat == quat:
    True
    >>> quat == quat.unit()
    False

or represent equivalent rotations:

    >>> quat.equal_rotation(quat.unit)
    True

Finally, the Quaternion class can efficiently handle an array of quaternions 
instead of just a single quaternion:

    >>> q = [[0.0, 0.0, 1.0, 1.0], [0.0, 1.0, 0.0, 1.0], [1.0, 0.0, 0.0, 1.0]]
    >>> quat = Quaternion(q)
    >>> quat
    <Quaternion (x, y, z, w)  [[0., 0., 1., 1.],
                               [0., 1., 0., 1.],
                               [1., 0., 0., 1.]] >

All attributes and methods work appropriately, for example:

    >>> quat * quat2
    Quaternion (x, y, z, w)  [[-1.,  1.,  1.,  1.],
                              [ 0.,  2.,  0.,  0.],
                              [ 1.,  1.,  1.,  1.]] >


The Spacecraft Frame
====================
There are several scenarios where it is necessary to know the position of the
spacecraft in orbit and its orientation.  For example, the spacecraft position
information is useful for determining if a source is occulted by the Earth at a
particular time, and the orientation is important for transforming between the
spacecraft and celestial frames. 

For this purpose, the |SpacecraftFrame| class defines the spacecraft position,
velocity, and a spacecraft rotation coordinate frame. This class inherits from 
Astropy's |BaseCoordinateFrame|_ and takes advantage of Astropy's 
transformation graph to enable transformation from a spacecraft frame to any 
celestial frame included in Astropy and vice versa.

Examples
--------
A |SpacecraftFrame| can handle either a defined position/velocity/orientation at
a specified time or an array of positions/velocities/orientations at a 
corresponding array of times.  To initialize a SpacecraftFrame, the inputs are
expected to be formatted as the following:

    *  ``obstime`` should be an Astropy Time object containing the time(s) at 
       which the position/velocity/orientation are determined;
    *  ``obsgeoloc`` should be an Astropy CartesianRepresentation representing 
       the spacecraft position in Earth-centered Cartesian Inertial coordinates;
    * ``obsgeovel`` should be an Astropy CartesianRepresentation representing
      the spacecraft velocity in Earth-centered Cartesian Inertial coordinates;
    * ``quaternion`` should be a |Quaternion| object.


    >>> from astropy.time import Time
    >>> from astropy.coordinates.representation import CartesianRepresentation
    >>> from gdt.core.coords import SpacecraftFrame, Quaternion
    >>>
    >>> times = Time(['2022-07-28 00:00:10', 
    >>>               '2022-07-28 00:00:11', 
    >>>               '2022-07-28 00:00:12'], format='iso')
    >>>
    >>> pos = CartesianRepresentation([(-5825973. ,  -5828639. , -5831298. ), 
    >>>                                ( 2209440. ,   2202317.2,  2195191.5), 
    >>>                                (-2965079. ,  -2965131.2, -2965179.5)],
    >>>                                unit='m')
    >>>
    >>> vel = CartesianRepresentation([(-2669.521 , -2662.4517,  -2655.377 ), 
    >>>                                (-7121.62  , -7124.296 ,  -7126.9634), 
    >>>                                (-53.883694, -50.27747 ,  -46.669613)],
    >>>                                unit='m/s')                                       
    >>>
    >>> quat = Quaternion([(0.41802019,  0.88160129, -0.01699675, -0.21851635),
    >>>                    (0.41793188,  0.88157926, -0.01744112, -0.21873903),
    >>>                    (0.41784446,  0.8815565 , -0.01788574, -0.21896175)])
    >>>
    >>> sc_frame = SpacecraftFrame(obstime=times, obsgeoloc=pos, obsgeovel=vel,
    >>>                            quaternion=quat)
    >>> sc_frame
    <SpacecraftFrame: 3 frames;
     obstime=[2022-07-28 00:00:10.000, ...]
     obsgeoloc=[(-5825973., 2209440., -2965079.) m, ...]
     obsgeovel=[(-2669.521, -7121.62, -53.883694) m / s, ...]
     quaternion=[(x, y, z, w) [ 0.41802019,  0.88160129, -0.01699675, -0.21851635], ...]>


The class is flexible in that if you don't have information for one of the 
inputs, then you can still instantiate an object with the other inputs. For 
example, if we don't have spacecraft velocity  information, we can initialize 
with only the position and orientation information.

SpacecraftFrame has a whole host of attributes that are inherited from Astropy's
BaseCoordinateFrame, so we will highlight the specialized attributes that the
SpacecraftFrame provides.  You can retrieve the location of the spacecraft in
terms of Earth latitude and longitude (and altitude):

    >>> sc_frame.earth_location.lat
    <Latitude [-25.7028489 , -25.70339896, -25.70391337] deg>
    
    >>> sc_frame.earth_location.lon
    <Longitude [-146.18754345, -146.12171964, -146.05589295] deg>
    
    >>> sc_frame.earth_location.height
    <Quantity [526242.35856409, 526239.37987702, 526236.25971664] m>

You can return the location of the spacecraft in the GCRS frame, which can then
be transformed to any other standard frame:

    >>> sc_frame.sc_gcrs
    <GCRS Coordinate (obstime=['2022-07-28 00:00:10.000' '2022-07-28 00:00:11.000'
     '2022-07-28 00:00:12.000'], obsgeoloc=(0., 0., 0.) m, obsgeovel=(0., 0., 0.) m / s): (ra, dec, distance) in (deg, deg, m)
        [(159.23128611, -25.44838127, 6900382.59711518),
         (159.30122099, -25.44887368, 6900379.45876531),
         (159.37115869, -25.44933074, 6900376.18928827)]>

Importantly, you can determine the location of the Earth geocenter as observed 
from the spacecraft and the apparent angular radius of the Earth:

    >>> sc_frame.geocenter
    <GCRS Coordinate (obstime=['2022-07-28 00:00:10.000' '2022-07-28 00:00:11.000'
     '2022-07-28 00:00:12.000'], obsgeoloc=(0., 0., 0.) m, 
     obsgeovel=(0., 0., 0.) m / s): (ra, dec, distance) in (deg, deg, m)
        [(339.23128611, 25.44838127, 6900382.59711518),
         (339.30122099, 25.44887368, 6900379.45876531),
         (339.37115869, 25.44933074, 6900376.18928827)]>
    
    >>> sc_frame.earth_angular_radius
    <Quantity [67.48524594, 67.48530557, 67.48536804] deg>
    
Finally, you can determine if the sun is visible, which utilizes the geocenter
and Earth angular radius information:

    >>> sc_frame.sun_visible
    array([ True,  True,  True])

In general, any location can be determined visible or not by the following:

    >>> from astropy.coordinates import SkyCoord
    >>> # RA=10.0, Dec=-30.0
    >>> coord = SkyCoord(10.0, -30.0, frame='icrs', unit='deg')
    >>> sc_frame.location_visible(coord)
    array([False, False, False])


We can also transform any celestial coordinate into our spacecraft frame:

    >>> sc_coord = coord.transform_to(sc_frame)
  
The coordinates of the SpacecraftFrame are azimuth (``az``) and elevation 
(``el``), and you can directly access the values through the SkyCoord object:

    >>> sc_coord.az
    <Longitude [123.49973775, 123.49278227, 123.48572411] deg>
    >>> sc_coord.el
    <Latitude [7.71285572, 7.66376691, 7.61466238] deg>
  
We can also do the inverse transformation: define a coordinate in the spacecraft
frame and transform it into the celestial frame.  If we wish to transform a 
single coordinate for all of the frames, then we must replicate that coordinate 
by the number of observation times we have in the SpacecraftFrame object:

    >>> # replicate the Azimuth and Elevation 3 times
    >>> coord = SkyCoord([37.]*3, [-5]*3, frame=sc_frame, unit='deg')
    >>> # transform to the ICRS frame
    >>> coord.icrs
    <SkyCoord (ICRS): (ra, dec) in deg
        [(88.18611282, 14.27154045), (88.19362827, 14.23809367),
         (88.20104308, 14.20460444)]>
  
    >>> # transform to the galactic frame
    >>> coord.galactic
    <SkyCoord (Galactic): (l, b) in deg
        [(193.4677391 , -6.08914974), (193.50052896, -6.0995751 ),
         (193.53330888, -6.11010218)]>

Any individual spacecraft frame from the object can be retrieved by index:

    >>> sc_frame[1]
    <SpacecraftFrame: 1 frames;
     obstime=[2022-07-28 00:00:11.000]
     obsgeoloc=[(-5828639., 2202317.2, -2965131.2) m]
     obsgeovel=[(-2662.4517, -7124.296, -50.27747) m / s]
     quaternion=[(x, y, z, w) [ 0.41793188,  0.88157926, -0.01744112, -0.21873903]]>

It carries all the same attributes and methods as the original object. For a 
single frame, you can transform multiple points from that frame into a different
frame:

    >>> coord = SkyCoord([37.0, 47.0, 57.0], [-5.0, -10.0, -15.0], 
                         frame=sc_frame[1], unit='deg')
    >>> coord.icrs
    <SkyCoord (ICRS): (ra, dec) in deg
        [(88.19362827, 14.23809367), (76.7368639 , 14.61547344),
         (65.40038798, 14.83838792)]>

Finally, if you have a SpacecraftFrame that contains multiple individual frames,
like in our example, you can interpolate and retrieve a frame at a given time:

    >>> requested_time = Time('2022-07-28 00:00:10.5', format='iso')
    >>> sc_frame.at(requested_time)
    <SpacecraftFrame: 1 frames;
     obstime=[2022-07-28 00:00:10.500]
     obsgeoloc=[(-5827306., 2205878.6, -2965105.1) m]
     obsgeovel=[(-2665.98635, -7122.958, -52.080582) m / s]
     quaternion=[(x, y, z, w) [ 0.41797605,  0.8815903 , -0.01721894, -0.2186277 ]]>


Designing a Class to Read From File
===================================
Most missions provide auxiliary data such as the spacecraft position and 
orientation along with other state information (whether the instruments were
observing, if the sun is visible, etc.).  It is then useful to be able to read
from a file and create the SpacecraftFrame from the data in the file. Below is 
an example of how you can design such a class to read from a FITS file.


    >>> from gdt.core.coords.spacecraft import SpacecraftFrameModelMixin, SpacecraftStatesModelMixin
    >>> from gdt.core.file import FitsFileContextManager
    >>>
    >>> class MyPosHistFile(SpacecraftFrameModelMixin, SpacecraftStatesModelMixin, FitsFileContextManager):
    >>> 
    >>>     def get_spacecraft_frame(self):
    >>>         # Format data in the file and create a SpacecraftFrame here.
    >>>
    >>>     def set_spacecraft_frame(self, spacecraft_frame):
    >>>         # Accept a SpacecraftFrame object so that the frame can be set
    >>>         # for e.g. writing the information to file.  This is optional.
    >>>
    >>>     def get_spacecraft_states(self):
    >>>         # Format data in the file and return an Astropy TimeSeries
    >>>
    >>>     def set_spacecraft_states(self, series):
    >>>         # Accept an Astropy TimeSerioes object so that the frame can be set
    >>>         # be set for e.g. writing the information to file.  
    >>>         # This is optional.

In this example, the |FitsFileContextManager| handles the FITS file access and
has some convenience functions (see :ref:`The File Module<core-file>` for 
details). The |SpacecraftFrameModelMixin| and the |SpacecraftStatesModelMixin| 
are abstract base classes defining the getter and setter for the frame and 
states, respectively.


Reference/API
=============

.. automodapi:: gdt.core.coords
   :inherited-members:
