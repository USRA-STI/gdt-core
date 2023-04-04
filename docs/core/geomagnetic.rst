.. _core-geomagnetic:
.. |SouthAtlanticAnomaly| replace:: :class:`~gdt.core.geomagnetic.SouthAtlanticAnomaly`


****************************************************
The Geomagnetic Module (:mod:`gdt.core.geomagnetic`)
****************************************************

Introduction
============
The Geomagnetic module primarily provides a base class for defining a 
South Atlantic Anomaly (SAA) polygon.  The SAA is a region of high radiation as
it is caused by the closest approach of the Van Allen belt to the Earth's 
surface.  Most gamma-ray detectors need to be deactivated whenever they pass
through the SAA, and the SAA boundary may be different depending on altitude and
observing energy range.  The |SouthAtlanticAnomaly| class is a container with
some convenience functions for an instrument's SAA boundary definition, 
although, in general, it can be used to define any other exclusionary region.

For Developers:
===============
To use the SouthAtlanticAnomaly class, we need to subclass it and 
define the latitude and longitude vertices of the boundary within the class.

    >>> from gdt.core.geomagnetic import SouthAtlanticAnomaly
    >>> 
    >>> class MySaaRegion(SouthAtlanticAnomaly):
    >>>     _latitude = [-30.000, -19.867, -9.733, 0.400, 2.000, 2.000, -1.000,
    >>>                  -6.155, -8.880, -14.220, -18.404, -30.000, -30.000]
    >>>     _longitude = [33.900, 12.398, -9.103, -30.605, -38.400, -45.000, -65.000,
    >>>                   -84.000, -89.200, -94.300, -94.300, -86.100, 33.900]

Here we set two mandatory class variables, ``_latitude`` and ``_longitude``, 
which represent the vertices of the polygon.  Then we can create an instance:

    >>> saa = MySaaRegion()
    >>> saa.latitude
    array([-30.   , -19.867,  -9.733,   0.4  ,   2.   ,   2.   ,  -1.   ,
            -6.155,  -8.88 , -14.22 , -18.404, -30.   , -30.   ])
    >>> saa.longitude
    array([ 33.9  ,  12.398,  -9.103, -30.605, -38.4  , -45.   , -65.   ,
           -84.   , -89.2  , -94.3  , -94.3  , -86.1  ,  33.9  ])

There are some convenience functions provided, such as determining if the 
region is a closed polygon:

    >>> saa.is_closed()
    True

and if the polygon is convex:

    >>> saa.is_convex()
    True

Reference/API
=============

.. automodapi:: gdt.core.geomagnetic
   :inherited-members:

