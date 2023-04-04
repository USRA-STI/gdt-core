.. _plot-plot:
.. |GdtPlot| replace:: :class:`~gdt.core.plot.plot.GdtPlot`
.. |PlotElement| replace:: :class:`~gdt.core.plot.plot.PlotElement`
.. |lib| replace:: :ref:`lib<plot-lib>`


*************************************************
Plot Element Classes (:mod:`~gdt.core.plot.plot`)
*************************************************

In the ``gdt.core.plot.plot`` module, there are two important base classes to be
aware of.  

*  |GdtPlot| is the base class for all plots in the GDT.  It 
   contains the basic functionality to instantiate plots and provide ways to 
   dynamically update them.
*  |PlotElement| is the base class representing a single plot element, such as
   a lightcurve histogram or a polygon plotted on the sky. These are generally
   wrappers around the lower-level |lib| functions with additional functions to 
   enable dynamic updating.

For the PlotElement base class, the basic plot attributes provided are the 
``alpha``, ``color``, ``visible``, and ``zorder``.  These attributes can be
updated dynamically.  Derived classes can add other attributes.  

An example of how to subclass PlotElement with an extra ``linestyle`` attribute:

    >>> from gdt.core.plot.plot import PlotElement
    >>> class MyPlotElement(PlotElement):
    >>>     def __init__(self, data, ax, color='C0', alpha=1.0, **kwargs):
    >>>         # data is some data that we want to plot
    >>>         # ax is the matplotlib axis we are plotting to
    >>>
    >>>         super().__init__(color=color, alpha=alpha)
    >>>         self._kwargs = kwargs
    >>>
    >>>         # artists is some iterable of matplotlib plot artists
    >>>         artists =  _create(data, ax)
    >>>
    >>>         # matplotlib is not consistent in how artists are represented,
    >>>         # so we should sanitize (or flatten) to a simple list with 
    >>>         # this base class method
    >>>         self._artists = self._sanitize_artists(artists)
    >>>
    >>>     def _create(self, data, ax):
    >>>         artists = function_calling_matplotlib(data, color=self.color,
    >>>                                               alpha=self.alpha, **self.kwargs)
    >>>         return artists
    >>>
    >>>     @property
    >>>     def linestyle(self):
    >>>         return self.artists[0].get_linestyle()
    >>>
    >>>     @linestyle.setter
    >>>     def linestyle(self, val):
    >>>         for artist in self.artists:
    >>>             artist.set_linestyle(val)


In this example, we store other keyword arguments internally, and then the 
actual plotting is performed in ``_create()``, where we create a matplotlib 
plotting object (or some collection of objects).  This example assumes that
``linestyle`` should be the same for all matplotlib artists contained in
``MyPlotElement``.

Reference/API
=============

.. automodapi:: gdt.core.plot.plot
   :inherited-members:


