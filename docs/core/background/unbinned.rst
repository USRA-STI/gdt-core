.. _background_unbinned:
.. |BackgroundFitter| replace:: :class:`~gdt.core.background.fitter.BackgroundFitter`
.. |NaivePoisson| replace:: :class:`~gdt.core.background.unbinned.NaivePoisson`


******************************************************************************
Background Algorithms for Unbinned Data (:mod:`~gdt.core.background.unbinned`)
******************************************************************************

Introduction
============
The :mod:`gdt.core.background.unbinned` module provides algorithms for background
fitting of unbinned (event) data.  Typically the background is fit in the time 
domain for each channel of multi-channel data. Currently an algorithm for 
estimating the background via a naive Poisson Maximum Likelihood method is 
provided.

The Naive Poisson Maximum Likelihood Algorithm
==============================================
This algorithm assumes the presence of time-varying Poisson process, and aims
to provide a MLE of the Poisson rate.  Since the data could contain bright 
sources, this is a naive approach to estimating the background and will fail
if bright sources are present.  The premise of the algorithm lies on the fact 
that the data can be binned at some minimum resolution such that each "bin"
contains either zero or one count.  The MLE Poisson solution can then be 
easily formulated by considering the product of zero-count and one-count 
Poisson likelihoods.  This formulation works for a constant Poisson process.
For a time-varying process, a sliding time-window concept is introduced such 
that within the time-window, the Poisson rate is assumed constant, and the 
Poisson MLE is made at the center of the window.  In a sense, this algorithm is
similar to "smoothing" the data as the window is moved through the data. 
Because the rate is estimated at the center of the window, when the window is
at either edges of the considered data, the window width is truncated, which 
effectively results in an increased uncertainty in the Poisson rate.

The |NaivePoisson| class contains an option between two versions of the
algorithm: a fast approximation, and a slower exact implementation.  The exact
version uses precisely a fixed window width as specified by the user (except
when truncated at the edges of the data), and counts the number of zero- and
one-count "bins".  The fast approximation allows the user-specified window
width to change by calculating the number of events contained within the first
window and keeping the number of events in the window constant instead of the 
window width.  For slowly varying data, this is a good approximation and runs
much faster than the exact implementation; however if the data rate changes 
considerably, the window width in this approximation will also change 
considerably.

Examples
--------
Although the algorithm is designed to fit multiple energy channels simultaneously,
for the sake of simplicity, we will look at an example using a single channel of
data.

    >>> import numpy as np
    >>> # generate events at an average rate of 1 per second for 1000 s.
    >>> times = np.random.exponential(1.0, size=1000).cumsum()

we have generated some random event data consistent with a constant rate of 
1 count/second. Now we can initialize the algorithm:

    >>> from gdt.core.background.unbinned import NaivePoisson
    >>> naivep = NaivePoisson([times])
    
Because we have a single channel of data, our array of times must be contained
in a single-element list.  Once initialized, we can do the fit.  There are two
algorithm-specific parameters we can set: the window width, and if we're using
the fast approximation or exact implementation:

    >>> naivep.fit(window_width=10.0, fast=False)

You can then interpolate the fit at any arbitrary point (within the data that
was fit):

    >>> interp_edges = np.linspace(10.0, 20.0, 6)
    >>> rate_interp, uncert_interp = naivep.interpolate(interp_edges[:-1], interp_edges[1:])
    >>> rate_interp
    array([[1.02945792],
           [0.82706984],
           [0.94514904],
           [1.04939443],
           [1.08009983]])
    >>> uncert_interp
    array([[0.32085167],
           [0.28758822],
           [0.30743276],
           [0.32394358],
           [0.32864872]])

You can refit the data with different parameters.  For example, here we use the
fast approximation:

    >>> naivep.fit(window_width=10.0, fast=True)
    >>> rate_interp, uncert_interp = naivep.interpolate(interp_edges[:-1], interp_edges[1:])
    >>> rate_interp
    array([[0.98921635],
           [0.97645178],
           [0.96820927],
           [0.99056896],
           [1.20206283]])
    >>> uncert_interp
    array([[0.2743579 ],
           [0.27067192],
           [0.26852942],
           [0.27461663],
           [0.33336566]])

.. _background_unbinned_plugins:
Designing Plug-ins
==================
You can design your own custom algorithm to be used to fit unbinned data.  In 
order for the algorithm to work with the fitter interface, the custom algorithm
must be a class and satisfy the following:

*  ``__init__()`` must take the following as an argument:
  #.  A list, of length ``num_chans``, where each item is a numpy.ndarray 
      of event times.
*  A ``fit()`` method that takes no arguments but may have keywords to specify
   algorithm parameters.
*  An ``interpolate()`` method that must take the following arguments:

   #.  A 1D array of time bin start times, shape (``num_times``,);
   #.  A 1D array of time bin end times, shape (``num_times``,).
*  The ``interpolate()`` method must return the following:

   #.  A 2D rates array, shape (``num_times``, ``num_chans``);
   #.  A 2D rate uncertainty array, shape (``num_times``, ``num_chans``).

Additionally, the class can provide the following public attributes to be 
exposed by the |BackgroundFitter| interface:

*  ``dof``: The degrees of freedom of the fits, array shape (``num_chans``,)
*  ``statistic``: The fit statistic for each fit, array shape (``num_chans``,)
*  ``statistic_name``: A string of the fit statistic used

To illustrate the above requirements, below is an example template:

    >>> class MyUnbinnedBackgroundAlgorithm():
    >>>
    >>>     def __init__(self, times):
    >>>         # Do some stuff here, like store the inputs for use by the algorithm
    >>>
    >>>     @property
    >>>     def dof(self):
    >>>         # Return the degrees of freedom here. This is optional.
    >>>
    >>>     @property
    >>>     def statistic(self):
    >>>         # Return the fit statistic for each channel.  This is optional.
    >>>
    >>>     @property
    >>>     def statistic_name(self):
    >>>         # Return the name of the fit statistic.  This is optional.
    >>>
    >>>     def fit(self, fit_param1=0, fit_param2=42.0):
    >>>         # Parameters needed for the algorithm are passed by keywords.
    >>>         # The algorithm for doing the fit goes here.
    >>>
    >>>     def interpolate(self, tstart, tstop, interp_param=10.0):
    >>>         # tstart, tstop required and extra parameters for interpolation
    >>>         # can be passed as keywords.
    >>>
    >>>         # The interpolation algorithm goes here.
    >>>         return (rates, rate_uncertainty)


Reference/API
=============

.. automodapi:: gdt.core.background.unbinned
   :inherited-members:
   :no-inheritance-diagram:
