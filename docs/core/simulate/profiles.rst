.. _sim-profiles:
.. |tophat| replace:: :func:`~gdt.core.simulate.profiles.tophat`

*******************************************************************
Source and Background Profiles (:mod:`~gdt.core.simulate.profiles`)
*******************************************************************
This package contains some functions for source pulse shapes and 
background time history trends.  While these are not complex, and the user can
certainly provide their own functions, these are provided for convenience and
as examples.  These profiles are used for generating simulated time-varying 
event data.

The general requirement for the profiles is that they must accept an array of
times and parameter values defining the profile, and they must return the 
function evaluated at the times.  For example:

    >>> def my_profile(times, param1, param2):
    >>>     # define function that evaluates profile at ``times``
    >>>     return func

Here, ``times`` are the times at which to evaluate the profile, ``param1``
and ``param2`` are parameters that define the profile, and ``func`` is the 
profile evaluated at ``times``.

As an example of how a profile is used, we can import the |tophat| profile:

    >>> import numpy as np
    >>> from gdt.core.simulate.profiles import tophat
    >>> times = np.linspace(0.0, 10.0, 11)
    
    >>> # amplitude = 1.0, tstart=5.0, tstop=8.0
    >>> tophat_params = (1.0, 5.0, 8.0)
    >>> tophat(times, *tophat_params)
    array([0., 0., 0., 0., 0., 1., 1., 1., 1., 0., 0.])
    

Reference/API
=============

.. automodapi:: gdt.core.simulate.profiles
   :inherited-members:


