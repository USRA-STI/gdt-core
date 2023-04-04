.. _binning_unbinned:
.. |numpy.histogram2d| replace:: numpy.histogram2d
.. _numpy.histogram2d: https://numpy.org/doc/stable/reference/generated/numpy.histogram2d.html


************************************************************************
Binning Algorithms for Unbinned Data (:mod:`~gdt.core.binning.unbinned`)
************************************************************************

A collection of algorithms are provided to perform binning on unbinned
(event) data.  Each algorithm is packaged in a function with the same
general inputs and outputs.  

For Developers:
===============
Users can define their own binning algorithm to be
used to bin event data within the data tools by following these rules:

*  The function must take as an argument ``times``, An array of event arrival times
*  Any algorithm-specific parameters can be pass as additional arguments or keywords
*  The function must return ``edges``, the array of bin edges.
  
Following this design, here is an example function:

    >>> def my_binning_algorithm(times, my_param1, my_param2=None):
    >>>     # my_param1 and my_param2 are algorithm-specific parameters
    >>>     # define algorithm here
    >>>     return edges

Note that the functions do not perform the binning but instead return the 
edges defining the bins.  The Data Tools uses the event times and the bin edges
with |numpy.histogram2d|_ (to handle multi-channel data) to do the binning once 
the edges are determined by the algorithm.

Reference/API
=============

.. automodapi:: gdt.core.binning.unbinned
   :inherited-members:
   :no-inheritance-diagram:

