.. _binning_binned:

********************************************************************
Binning Algorithms for Binned Data (:mod:`~gdt.core.binning.binned`)
********************************************************************

A collection of algorithms are provided to perform binning/rebinning on binned
(histogram) data.  Each algorithm is packaged in a function with the same
general inputs and outputs.  

For Developers:
===============
Users can define their own binning algorithm to be
used to bin data within the data tools by following these rules:

*  The function must take as arguments:
  
  #. ``counts``: An array of counts in each bin, size `m`
  
  #.  ``exposure``: An array of the exposure of each bin, size m
  
  #.  ``old_edges``: An array of the current bin edges, size `m` + 1
  
*  Any algorithm-specific parameters can be pass as additional arguments or keywords

*  The function must return:
  
  #.  ``new_counts``: The array of counts in each of the new bins, size `n`
  
  #.  ``new_exposure``: The array of exposure for each of the new bins, size `n`
  
  #.  ``new_edges``: The array of new bin edges, size `n` + 1
  
Following this design, here is an example function:

    >>> def my_binning_algorithm(counts, exposure, old_edges, my_param1):
    >>>     # my_param1 is an algorithm-specific parameter
    >>>     # define algorithm here
    >>>     return (new_counts, new_exposure, new_edges)


Reference/API
=============

.. automodapi:: gdt.core.binning.binned
   :inherited-members:
   :no-inheritance-diagram:
