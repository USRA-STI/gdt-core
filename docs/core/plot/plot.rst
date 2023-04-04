.. _plot:
.. |lib| replace:: :ref:`lib<plot-lib>`
.. |plot| replace:: :ref:`plot<plot-plot>`

***************************************
The Plot Package (:mod:`gdt.core.plot`)
***************************************

The ``plot`` package contains generalized plot functions and classes for 
visualization in the GDT.  Generally, there are three levels of abstraction in
this package:

#. Low-level plotting functions contained in |lib|;
#. Mid-level plot element classes that provide an interface to the the low-level 
   functions and allow for dynamic update of plots, contained in |plot|
#. High-level plot classes that group together a number of mid-level plot
   elements into a figure.

The package is divided into the following modules:

.. toctree::
   :maxdepth: 1

   drm
   earthplot
   lightcurve
   model
   plot_plot
   sky
   spectrum
   lib