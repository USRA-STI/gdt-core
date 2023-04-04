.. _core-tte:
.. |PhotonList| replace:: :class:`~gdt.core.tte.PhotonList`
.. |Phaii| replace:: :class:`~gdt.core.phaii.Phaii`
.. |Pha| replace:: :class:`~gdt.core.pha.Pha`
.. |EnergyBins| replace:: :class:`~gdt.core.data_primitives.EnergyBins`
.. |EventList| replace:: :class:`~gdt.core.data_primitives.EventList`
.. |Gti| replace:: :class:`~gdt.core.data_primitives.Gti`
.. |Ebounds| replace:: :class:`~gdt.core.data_primitives.Ebounds`
.. |PhotonList.open()| replace:: :meth:`~gdt.core.tte.PhotonList.open`
.. |PhotonList.write()| replace:: :meth:`~gdt.core.tte.PhotonList.write`
.. |PhotonList._build_hdulist()| replace:: :meth:`~gdt.core.tte.PhotonList._build_hdulist`
.. |PhotonList._build_headers()| replace:: :meth:`~gdt.core.tte.PhotonList._build_headers`
.. |FileHeaders| replace:: :class:`~gdt.core.headers.FileHeaders`
.. |core-headers| replace:: :ref:`Data File Headers<core-headers>`

**************************************************************
Photon List and Time-Tagged Event Files  (:mod:`gdt.core.tte`)
**************************************************************

Introduction
============
While the :ref:`PHAII Documentation<core-phaii>` describe a *binned* time-series 
of spectra, Photon List or Time-Tagged Event (TTE) data represent an *unbinned* 
time series.  In practice, this type of data is only semi-unbinned in the sense 
that while the time axis of the data are unbinned (hence the time tags), the 
energy axis is necessarily binned. This is because the energy resolution and 
response of an instrument does not allow the relative precision that is easily 
achievable for the time axis.  Therefore, the model for Photon List/TTE data is 
a series of photons/counts, each with at least two attributes: a time and an 
energy channel.

The PhotonList Class
====================
The |PhotonList| class provides a way to construct, write out, and read photon 
list and TTE files. There are also some functions provided to operate on this
data. Similar to the |Phaii| class, the PhotonList base class does not have the 
ability to read and write files.  PhotonList is expected to be subclassed to
define the reading and writing portion of the file, however, you can still
programmatically create a PhotonList object.  In the following examples, we will 
first create a PhotonList object with some data, and then we will walk through 
how to subclass PhotonList to create a class that can read/write files.

Examples
--------
First, we will show how to construct a |PhotonList| object from data 
count spectra. The data within a PhotonList object is an |EventList|, which is a
list of times and associated energy channel numbers, so we will define that. 
Additionally, we can define a |Gti|, which is one or multiple Good Time 
Intervals over which the data can be used for science.

  >>> import numpy as np
  >>> from gdt.core.data_primitives import EventList, Ebounds, Gti
  >>> # simulated Poisson rate of 1 count/sec
  >>> times = np.random.exponential(1.0, size=100).cumsum()
  >>> # random channel numbers
  >>> channels = np.random.randint(0, 6, size=100)
  >>> # channel-to-energy mapping
  >>> ebounds = Ebounds.from_bounds([10.0, 20.0, 40.0, 80.0, 160.0, 320.0], 
  >>>                               [20.0, 40.0, 80.0, 160.0, 320.0, 640.0])
  >>> data = EventList(times, channels, ebounds=ebounds)
  
  >>> # construct the good time interval(s)
  >>> gti = Gti.from_list([(0.0000, 100.0)])

  >>> # create the PhotonList object
  >>> from gdt.core.tte import PhotonList
  >>> tte = PhotonList.from_data(data, gti=gti, trigger_time=356223561.,
                                 event_deadtime=0.001, overflow_deadtime=0.1)
  >>> tte
  <PhotonList: 
   trigger time: 356223561.0;
   time range (1.2671749481077967, 94.56816998354927);
   energy range (10.0, 640.0)>

When the PhotonList object is created, note that we specified ``event_deadtime``
and ``overflow_deadtime``, which are the per-event deadtime imposed for the 
energy channels and overflow channel (assumed to be the last energy channel), 
respectively. Defining these variables are necessary to accurately determine
the exposure for a variety of functions.  If these variables are not set, then
it is assumed that there is no deadtime to be applied.

Additionally, setting the GTI is not required, and omitting the GTI will cause a 
default GTI to be created with a range spanning the full time range of the data.  
Also note that we set a trigger time (also optional), that is used as a time 
reference for the data.

Now that we have created our PhotonList object, there are several attributes 
that are available to us.  We can directly access the data and GTI we created 
the object with:

  >>> # the EventList data
  >>> tte.data
  <EventList: 100 events;
   time range (1.2671749481077967, 94.56816998354927);
   channel range (0, 5)>
  >>> # the GTI
  >>> tte.gti
  <Gti: 1 intervals; range (0.0, 100.0)>
  >>> # the Ebounds
  >>> tte.ebounds
  <Ebounds: 6 intervals; range (10.0, 640.0)>

There are other attributes that are exposed:

  >>> # energy range
  >>> tte.energy_range
  (10.0, 640.0)
  >>> # number of energy channels
  >>> tte.num_chans
  6
  >>> # time range covered by data (time of first and last events)
  >>> tte.time_range
  (1.2671749481077967, 94.56816998354927)
  >>> # trigger time (if available)
  >>> tte.trigtime
  356223561.0

In addition to these attributes, we can retrieve the exposure for the full data
contained within the PhotonList object or for a time segment of the data 
contained:

  >>> # full exposure
  >>> tte.get_exposure()
  91.518

  >>> # get total exposure for two segments of data
  >>> tte.get_exposure(time_ranges=[(0.0, 10.0), (20.0, 30.0)])
  19.786

You can slice the PhotonList object in time or energy.  You can specify one or
multiple ranges to slice over:

  >>> sliced_energy = tte.slice_energy([(25.0, 35.0), (550.0, 600.0)])
  >>> sliced_energy.data
  <EventList: 35 events;
   time range (3.3760599925164128, 94.56816998354927);
   channel range (1, 5)>

  >>> sliced_time = tte.slice_time([(0.0, 10.0), (20.0, 30.0)])
  >>> sliced_time.data
  <EventList: 16 events;
   time range (1.2671749481077967, 29.060934691909935);
   channel range (0, 5)>

In both cases, we sliced to two disjoint ranges, defined as a list of tuples.  
Note that the resulting sliced energy edges are dependent on the energy edges 
of the original object, since they cannot be changed.


Similar to the PHA and PHAII data, a PhotonList can be rebinned in energy.  Here
we rebin the energy channels by a factor of 2:

  >>> from gdt.core.binning.binned import combine_by_factor
  >>> rebinned_energy = tte.rebin_energy(combine_by_factor, 2)
  >>> rebinned_energy.data
  <EventList: 100 events;
   time range (1.2671749481077967, 94.56816998354927);
   channel range (0, 2)>

Since the time axis of the data is unbinned, the approach to binning it is a
little different than that for PHAII data.  In fact, what you can do is directly
convert the PhotonList object to a Phaii object given a binning function:

  >>> from gdt.core.binning.unbinned import bin_by_time
  >>> from gdt.core.phaii import Phaii
  >>> phaii = tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii)
  >>> phaii
  <Phaii: 
   trigger time: 356223561.0;
   time range (1.2671749481077967, 95.26717494810778);
   energy range (10.0, 640.0)>

There are a couple things to note here.  In this example, we binned the data in
time (using ``bin_by_time``) to 1.0 s resolution.  The other thing to note is 
that we set the keyword argument ``phaii_class`` to the base Phaii class.  If
you are subclassing both PhotonList and Phaii for a particular type of data
from an instrument, you can assign your subclassed Phaii instead 
(e.g. ``phaii_class=MyPhaii``).  This will apply any specialized header and HDU
construction to your output Phaii object. Finally, you can choose to only
convert a segment of the PhotonList to a Phaii:

  >>> # convert a time range
  >>> tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, time_range=(0.0, 10.0))
  <Phaii: 
   trigger time: 356223561.0;
   time range (1.2671749481077967, 10.267174948107796);
   energy range (10.0, 640.0)>
  
  >>> # convert an energy range
  >>> tte.to_phaii(bin_by_time, 1.0, phaii_class=Phaii, energy_range=(50, 300))
  <Phaii: 
   trigger time: 356223561.0;
   time range (1.2671749481077967, 92.26717494810778);
   energy range (40.0, 320.0)>

Similar to Phaii data, a PhotonList can be integrated over time to produce a 
count spectrum, represented as an |EnergyBins| object:

  >>> # integrate over the full time range
  >>> spectrum = tte.to_spectrum()
  <EnergyBins: 6 bins;
   range (10.0, 640.0);
   1 contiguous segments>

  >>> # integrate over a subset of the full time range
  >>> spectrum = tte.to_spectrum(time_range=(0.0, 10.0))

We can even directly create a |Pha| object, which integrates over a time range
(or multiple time ranges) and produces a fully qualified object that can then
be written to a FITS file.

  >>> pha = tte.to_pha(time_ranges=[(0.0, 10.0), (20.0, 30.0)])
  >>> pha
  <Pha: 
   trigger time: 356223561.0;
   time range (1.2671749481077967, 29.060934691909935);
   energy range (10.0, 640.0)>
  >>> pha.gti
  <Gti: 2 intervals; range (1.2671749481077967, 29.060934691909935)>

Furthermore, you can create a Pha object with only a subset of the energy range,
and it will automatically mask out the channels you are not using:

  >>> pha = tte.to_pha(energy_range=(50.0, 300.0))
  <Pha: 
   trigger time: 356223561.0;
   time range (1.2671749481077967, 94.56816998354927);
   energy range (10.0, 640.0)>
  >>> pha.channel_mask
  array([False, False,  True,  True,  True, False])
  >>> pha.valid_channels
  array([2, 3, 4])

Finally, you can merge multiple PhotonList objects together into a single 
object.  A word of caution: real data contained within a PhotonList may be
quite large, so make sure you have sufficient memory to perform such a merge.

  >>> tte1 = tte.slice_time((0.0, 10.0))
  >>> tte2 = tte.slice_time((20.0, 30.0))
  >>> tte_merged = PhotonList.merge([tte1, tte2])
  >>> tte_merged
  <PhotonList: 
   trigger time: 356223561.0;
   time range (1.2671749481077967, 29.060934691909935);
   energy range (10.0, 640.0)>
  >>> phaii_merged.gti
  <Gti: 2 intervals; range (1.2671749481077967, 29.060934691909935)>


Subclassing
-----------
To read and write PhotonList/TTE FITS files, the |PhotonList| class needs to be 
subclassed. This is because the format and metadata of the files can be 
different from mission to mission.  When subclassing PhotonList to read a 
file, the |PhotonList.open()| method needs to be defined.  To write out a file, 
the private method |PhotonList._build_hdulist()| needs to be defined, which 
defines the data structure for each extension of the FITS file. Adding header 
information/metadata is not required, however if you do want the header 
information to be tracked when reading in a file and written out when writing a
file to disk, you will need to create the header definitions as explained in
|core-headers| and also define the private method |PhotonList._build_headers()|.

To illustrate further, below is a sketch of how the PhotonList class should be 
subclassed in the example ``MyPhotonList``:

    >>> import astropy.io.fits as fits
    >>> from gdt.core.data_primitives import Ebounds, EventList, Gti
    >>> from gdt.core.tte import PhotonList
    >>>
    >>> class MyPhotonList(PhotonList):
    >>>     """An example to read and write PhotonList/TTE files for xxx instrument"""
    >>>     @classmethod
    >>>     def open(cls, file_path, **kwargs):
    >>>         with super().open(file_path, **kwargs) as obj:
    >>>
    >>>             # an example of how to set the headers
    >>>             hdrs = [hdu.header for hdu in obj.hdulist]
    >>>             headers = MyFileHeaders.from_headers(hdrs)
    >>>
    >>>             # an example of how to set the ebounds
    >>>             ebounds = Ebounds.from_list(emin, emax)
    >>>
    >>>             # an example of how to set the data
    >>>             data = EventList(times, channels, ebounds=ebounds)
    >>>
    >>>             # an example of how to set the GTI
    >>>             gti = Gti.from_bounds(gti_start, gti_end)
    >>>
    >>>             return cls.from_data(data, gti=gti, trigger_time=trigger_time,
    >>>                                  filename=obj.filename,czheaders=headers)
    >>>
    >>>     def _build_hdulist(self):
    >>>         # this is where we build the HDU list (header/data for each extension)
    >>>         hdulist = fits.HDUList()
    >>>
    >>>         # some code to create PRIMARY HDU
    >>>         # ...
    >>>         hdulist.append(primary_hdu)
    >>>
    >>>         # code to create other HDUs and append to hdulist
    >>>         # ...
    >>>
    >>>         return hdulist
    >>>
    >>>     def _build_headers(self, trigtime, tstart, tstop, num_chans):
    >>>         # build the header based on these inputs
    >>>         headers = self.headers.copy()
    >>>         # update header information here
    >>>         # ...
    >>>
    >>>         return headers


To create a PhotonList object from a FITS file, the ``open()`` method should, 
at a minimum, be able to construct an |EventList| object containing the 
data.  Additionally, you can construct a |Gti|, and |FileHeaders|. If the data
has an associated trigger or reference time, you can track that as well.

The example creation of ``headers`` takes in a list of the headers from
each extension and assumes you have created a class (in this case 
``MyFileHeaders``) that will read in the header information.

To write the data to disk, the |PhotonList._build_hdulist()| defines the list of 
FITS HDUs, and is called by the |PhotonList.write()| method.  The 
|PhotonList._build_headers()| method is called whenever operations are performed 
on the object, like rebinning or slicing to propagate the header information.


Reference/API
=============

.. automodapi:: gdt.core.tte
   :inherited-members:

Special Methods
===============
.. automethod:: gdt.core.tte.PhotonList._build_hdulist
.. automethod:: gdt.core.tte.PhotonList._build_headers

