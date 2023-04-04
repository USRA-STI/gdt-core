.. _core-time:
.. |Rsp| replace:: :class:`~gdt.core.response.Rsp`
.. |Rsp2| replace:: :class:`~gdt.core.response.Rsp2`

*************************************************
Time Epochs and Utilities  (:mod:`gdt.core.time`)
*************************************************
The GDT relies heavily on Astropy's infrastructure for time definitions and 
conversions.

For Developers:
===============

Custom Time Epochs
------------------
Missions often define a *Mission Elapsed Time* (MET), which is essentially a 
time epoch that serves as a reference for all observations and data produced
during the lifetime of the mission.  METs can be defined in different ways, 
on different time scales, and with different reference points.  Astropy 
provides a way to define our own custom epoch so that we can convert between
our custom MET to any standard time scale and convention, or even to other 
custom METs.  Defining a custom epoch is not difficult, and we show, as an 
example, the definition for the Fermi mission epoch:

    >>> from astropy.time import TimeFromEpoch
    >>>
    >>> class FermiSecTime(TimeFromEpoch):
    >>>     """Represents the number of seconds elapsed since Jan 1, 2001 00:00:00 UTC
    >>>     including leap seconds.
    >>>     """
    >>>     
    >>>     name = 'fermi'
    >>>     unit = 1.0 / 86400.0  # in days (1 day == 86400 seconds)
    >>>     epoch_val = '2001-01-01 00:01:04.184' # offset of TT from UTC
    >>>     epoch_val2 = None # Do not need epoch_val2
    >>>     epoch_format = 'iso'
    >>>     epoch_scale = 'tt'  # Terrestrial Time
    
In this example, there a few things to note. First the class must inherit 
Astropy's ``TimeFromEpoch``, and the name of our epoch is captured in the 
variable ``name``.  Next, ``unit`` defines the units of time relative to a day.
For example, if we expect the unit of time to be 1 second, then ``unit`` is
defined as it is in the example above.  The other variables define the reference
point and scale for the custom epoch.

Examples
--------
The custom epoch can be used in the following way.  For example, let us create
a UTC time in an ISO format:

    >>> from astropy.time import Time
    >>> time = Time('2022-01-01 12:22:42', format='iso', scale='utc')

Now, the time can be converted to our custom epoch by using the name we gave 
epoch:

    >>> time.fermi
    662732567.0

That is the corresponding Fermi MET for the given UTC time.  We can go the other
way, by defining the Fermi MET in a time object:

    >>> time = Time(662732567.0, format='fermi')
    >>> time
    <Time object: scale='tt' format='fermi' value=662732567.0>

And we can convert it back to the original UTC scale and ISO format:
    
    >>> time.utc.iso
    '2022-01-01 12:22:42.000'

The time module contains some utilities, such as generating a range of times:

    >>> from gdt.core.time import time_range
    >>> # generate times from our MET input time, to +100 s in increment of 10 s
    >>> times = time_range(time, 100.0, 10)
    >>> times
    <Time object: scale='tt' format='fermi' value=[6.62732567e+08 6.62732577e+08 6.62732587e+08 6.62732597e+08
     6.62732607e+08 6.62732617e+08 6.62732627e+08 6.62732637e+08
     6.62732647e+08 6.62732657e+08]>
    >>> times.utc.iso
    array(['2022-01-01 12:22:42.000', '2022-01-01 12:22:52.000',
           '2022-01-01 12:23:02.000', '2022-01-01 12:23:12.000',
           '2022-01-01 12:23:22.000', '2022-01-01 12:23:32.000',
           '2022-01-01 12:23:42.000', '2022-01-01 12:23:52.000',
           '2022-01-01 12:24:02.000', '2022-01-01 12:24:12.000'], dtype='<U23')



Reference/API
=============

.. automodapi:: gdt.core.time
   :inherited-members:


