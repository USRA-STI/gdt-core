.. _background_binned:
.. |BackgroundFitter| replace:: :class:`~gdt.core.background.fitter.BackgroundFitter`


**************************************************************************
Background Algorithms for Binned Data (:mod:`~gdt.core.background.binned`)
**************************************************************************

Introduction
============
The :mod:`gdt.core.background.binned` module provides algorithms for background
fitting of binned data.  Typically the background is fit in the time domain
for each channel of multi-channel data, however, nothing precludes the design 
of an algorithm that fits in the energy domain. Currently an algorithm for 
fitting a polynomial to binned data is provided.

The Polynomial Fitting Algorithm
================================
The algorithm for fitting a polynomial to binned data is based on least-squares
minimization using chi-squared as a statistic.  The fit is performed in two 
passes to mitigate issues with small count statistics.  The first pass uses 
variances calculated from the observed data rates in the minimization. The
second pass uses variances calculated from the first-pass model rates.  The fit
is performed on multi-channel data, where each channel is fit independently with
the same order polynomial.

Examples
--------
Although the algorithm is designed to fit multiple energy channels simultaneously,
for the sake of simplicity, we will look at an example using a single channel of
data.

    >>> import numpy as np
    >>> # generate events at an average rate of 1 per second for 1000 s.
    >>> times = np.random.exponential(1.0, size=1000).cumsum()
    
    >>> from gdt.core.binning.unbinned import bin_by_time
    >>> # get the bin edges such that we have equally spaced 10-sec wide bins
    >>> edges = bin_by_time(times, 10.0)
    
    >>> counts, _ = np.histogram(times, bins=edges)

So far we have generated some random event data consistent with a constant
rate of 1 count/second and binned it into 10 s bins.  Because this is a single
channel of data, and Polynomial expects multiple channels, we just need to 
reshape the counts array.  Furthermore, let's assume our exposure is the same
as the bin width.

    >>> counts = counts.reshape(-1, 1)
    >>> exposure = edges[1:] - edges[:-1]
    
Now we can initialize the Polynomial object:

    >>> from gdt.core.background.binned import Polynomial
    >>> poly = Polynomial(counts, edges[:-1], edges[1:], exposure)

Once the object is initialized, we can do the fit and pass any algorithms 
parameters needed to fully specify the fit. In this example, we will pass the
polynomial order.

    >>> model, model_uncert = poly.fit(order=0)
    >>> print(model)
    [[1.02040816]
     [1.02040816]
     [1.02040816]
    ...
    >>> print(model_uncert)
    [[0.03019387]
     [0.03019387]
     [0.03019387]
    ...

The fit returns the model fit at each bin as well as the model uncertainty.  
Since this is based off of randomly generated data, your results might be 
slightly different.  You can also retrieve the fit statistic and 
degrees-of-freedom for each channel of the fit:

    >>> poly.statistic
    array([97.992])
    >>> poly.dof
    array([97.])

Of course, this is a trivial case, so maybe we want to fit a higher-order
polynomial and see if it provides a good fit.  You can simply redo the fit:

    >>> model, model_uncert = poly.fit(order=1)
    >>> poly.statistic
    array([97.86047063])
    >>> poly.dof
    array([96.])

Finally, we can interpolate the rate and rate uncertainty at any point:

    >>> interp_edges = np.linspace(10.0, 20.0, 6)
    >>> model_interp, uncert_interp = poly.interpolate(interp_edges[:-1], interp_edges[1:])
    >>> model_interp
    array([[1.00572787],
           [1.00578859],
           [1.00584931],
           [1.00591003],
           [1.00597075]])
    >>> uncert_interp
    array([[0.05940458],
           [0.05922084],
           [0.05903729],
           [0.05885395],
           [0.05867081]])


.. _background_binned_plugins:
Designing Plug-ins
==================
You can design your own custom algorithm to be used to fit binned data.  In 
order for the algorithm to work with the fitter interface, the custom algorithm
must be a class and satisfy the following:

*  ``__init__()`` must take the following as arguments:
  #.  A 2D array of counts with shape (``num_times``, ``num_chans``);
  #.  A 1D array of time bin start times, shape (``num_times``,); 
  #.  A 1D array of time bin end times, shape (``num_times``,);
  #.  A 1D array of exposures for each time bin, shape (``num_times``,).
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

    >>> class MyBinnedBackgroundAlgorithm():
    >>>
    >>>     def __init__(self, counts, tstart, tstop, exposure):
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

.. automodapi:: gdt.core.background.binned
   :inherited-members:
   :no-inheritance-diagram:


