.. _background_primitives:
.. |BackgroundRates| replace:: :class:`~gdt.core.background.primitives.BackgroundRates`
.. |BackgroundSpectrum| replace:: :class:`~gdt.core.background.primitives.BackgroundSpectrum`
.. |TimeEnergyBins| replace:: :class:`~gdt.core.data_primitives.TimeEnergyBins`
.. |EnergyBins| replace:: :class:`~gdt.core.data_primitives.EnergyBins`
.. |BAK| replace:: :class:`~gdt.core.pha.Bak`

**************************************************************************
Background Primitive Data Classes (:mod:`~gdt.core.background.primitives`)
**************************************************************************

Introduction
=============
There are two "primitive" data classes for background data: the 
|BackgroundRates| class, which contains background model rates and uncertainties
for a time history of multi-channel data, and the |BackgroundSpectrum| class, 
which contains the background model count spectrum.  The former largely inherits
from the |TimeEnergyBins| primitive, while the latter inherits from the 
|EnergyBins| primitive.

The BackgroundRates Class
==========================
The model for the BackgroundRates class is the same as for |TimeEnergyBins|
with two important distinctions:

#.  The initialization of BackgroundRates requires **rates** instead of 
    **counts**. This is because we are storing a model instead of observed counts.
#.  The initialization also requires a **rate uncertainty** because the 
    background model may have Gaussian (or other) uncertainties associated with it.
    
Those differences on initialization aside, BackgroundRates operates, and
has most of the same attributes and methods that |TimeEnergyBins| has.

Examples
--------
As an example, let's assume we have a segment of a background model.  We need
the modeled background rate and uncertainty in each channel and each time bin. 
Of course we also need to know the bin edges and exposures.

    >>> from gdt.core.background.primitives import BackgroundRates
    >>> # three energy channels and four time bins
    >>> rates = [[30.0, 50.0, 10.0], [35.0, 55.0, 15.0], 
    >>>          [40.0, 60.0, 20.0], [45.0, 65.0, 25.0]]
    >>> # 10% uncertainty
    >>> rate_uncert = [[3.0, 5.0, 1.0], [3.5, 5.5, 1.5], 
    >>>                [4.0, 6.0, 2.0], [4.5, 6.5, 2.5]]
    >>>
    >>> tstart = [0.0, 1.0, 2.0, 3.0]
    >>> tstop = [1.0, 2.0, 3.0, 4.0]
    >>> exposure = exposure = [1.0] * 4
    >>> emin = [10.0, 50.0, 150.0]
    >>> emax = [50.0, 150., 300.0]
    >>> back_rates = BackgroundRates(rates, rate_uncert, tstart, tstop, emin, emax, 
    >>>                              exposure=exposure)
    >>> back_rates
    <BackgroundRates: 4 time bins;
     time range (0.0, 4.0);
     1 time segments;
     3 energy bins;
     energy range (10.0, 300.0);
     1 energy segments>

In our simple example, we've assumed each rate has a 10% relative uncertainty,
and that we have no deadtime (exposure equals bin width).  You can access all
of the same attributes as you would with |TimeEnergyBins|.  For example

    >>> back_rates.time_centroids
    array([0.5, 1.5, 2.5, 3.5])
    
    >>> back_rates.rates_per_kev
    array([[0.75      , 0.5       , 0.06666667],
           [0.875     , 0.55      , 0.1       ],
           [1.        , 0.6       , 0.13333333],
           [1.125     , 0.65      , 0.16666667]])

Also similar to |TimeEnergyBins| you can slice by energy or time:

    >>> back_rates.slice_energy(20.0, 100.0)
    <BackgroundRates: 4 time bins;
     time range (0.0, 4.0);
     1 time segments;
     2 energy bins;
     energy range (10.0, 150.0);
     1 energy segments>
    
    >>> back_rates.slice_time(1.0, 3.0)
    <BackgroundRates: 2 time bins;
     time range (1.0, 3.0);
     1 time segments;
     3 energy bins;
     energy range (10.0, 300.0);
     1 energy segments>
     
You can merge multiple BackgroundRates objects in time:
    
    >>> back_slice1 = back_rates.slice_time(0.0, 1.0)
    >>> back_slice2 = back_rates.slice_time(3.0, 4.0)
    >>> BackgroundRates.merge_time([back_slice1, back_slice2])
    <BackgroundRates: 2 time bins;
     time range (0.0, 4.0);
     2 time segments;
     3 energy bins;
     energy range (10.0, 300.0);
     1 energy segments>

Note that we slice the object into two discontiguous segments and then merged 
them together.  You can even sum multiple BackgroundRates together along the 
time axis, like in a scenario where you want to add the backgrounds from 
multiple detectors together:

    >>> back_rates_double = BackgroundRates.sum_time([back_rates, back_rates])
    >>> back_rates_double.rate
    array([[ 60., 100.,  20.],
           [ 70., 110.,  30.],
           [ 80., 120.,  40.],
           [ 90., 130.,  50.]])
    
You can integrate a BackgroundRates object over energy to create a new, 
single-channel BackgroundRates object: 

    >>> # integrate over full energy range
    >>> back_rates.integrate_energy()
    <BackgroundRates: 4 time bins;
     time range (0.0, 4.0);
     1 time segments;
     1 energy bins;
     energy range (10.0, 300.0);
     1 energy segments>

    >>> # integrate over partial energy range
    >>> back_rates.integrate_energy(emin=50., emax=100.0)
    <BackgroundRates: 4 time bins;
     time range (0.0, 4.0);
     1 time segments;
     1 energy bins;
     energy range (50.0, 150.0);
     1 energy segments>
    
Finally, the BackgroundRates can be integrated over time to produce a 
|BackgroundSpectrum| object:

    >>> back_rates.integrate_time()
    <BackgroundSpectrum: 3 bins;
     range (10.0, 300.0);
     1 contiguous segments>

Or you can directly create a background |BAK| object containing metadata that
can the be written to disk:
    
    >>> back_rates.to_bak(time_range=[1.0, 2.0], trigger_time=12345.678,
                          filename='holla.bak')
    <Bak: holla.bak;
     trigger time: 12345.678;
     time range (1.0, 2.0);
     energy range (10.0, 300.0)>


The BackgroundSpectrum Class
============================
As with BackgroundRates, two distinctions separate the initialization of 
|BackgroundSpectrum| from |EnergyBins|:

#.  Instead of initializing with counts, BackgroundSpectrum is initialized
    with the background rate in each energy channl.
#.  BackgroundSpectrum must also be initialized with the rate uncertainty in
    each energy channel.

Other than those differences in initialization, most of the attributes and 
methods are shared with |EnergyBins|.

Examples
--------
Often BackgroundSpectrum isn't initialized by the user, but rather is derived
from BackgroundRates; however it can be initialized in the following way:

    >>> from gdt.core.background.primitives import BackgroundSpectrum
    >>> rates = [37.5, 57.5, 17.5]
    >>> rate_uncert = [1.896, 2.889, 0.919]
    >>>
    >>> emin = [10.0, 50.0, 150.0]
    >>> emax = [50.0, 150., 300.0]
    >>> exposure = 4.0
    >>> back_spec = BackgroundSpectrum(rates, rate_uncert, emin, emax, exposure)
    >>> back_spec
    <BackgroundSpectrum: 3 bins;
     range (10.0, 300.0);
     1 contiguous segments>

The BackgroundSpectrum object has a variety of attributes, including the counts
and count uncertainty in each channel, which are derived quantities:

    >>> back_spec.counts
    array([150., 230.,  70.])
    >>> back_spec.count_uncertainty
    array([ 7.584, 11.556,  3.676])
    
The BackgroundSpectrum object can be sliced in energy in the following way:

    >>> back_spec.slice(50.0, 100.0)
    <BackgroundSpectrum: 1 bins;
     range (50.0, 150.0);
     1 contiguous segments>
     
Multiple BackgroundSpectrum objects can be also be summed, as long as they have
the same energy range and number of channels:
    
    >>> back_spec_double = BackgroundSpectrum.sum([back_spec, back_spec])
    >>> back_spec_double.counts
    array([300., 460., 140.])


Reference/API
=============

.. automodapi:: gdt.core.background.primitives
   :inherited-members:


