.. _gdt-core:

Welcome to the Gamma-ray Data Tools Core Package!
=================================================

.. figure:: images/gdt_logo_big.png

The Gamma-ray Data Tools (GDT) is centralized toolkit for hard X-ray and 
gamma-ray astrophysics data analysis, with a focus on providing a uniform 
interface to the data provided by several different missions and instruments.

The GDT Core Package (``gdt-core``) contains the core components of the GDT that
can be utilized for various instruments and is a generalized version of the
`Fermi GBM Data Tools <https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs>`_.
Individual mission or instrument packages can be developed using the ``gdt-core`` 
and released under the ``gdt`` namespace (see ``gdt-fermi`` as an example).

The documentation linked below walks through all of sub-packages and modules
within the gdt-core, and developers should take special note of the 
"**For Developers**" sections that detail how to subclass or design plugin functions
or classes for use in new instrument packages.

.. rubric:: Citing

If you use the GDT Core package to develop your own mission our instrument
package, we would appreciate an appropriate acknowledgment. For publications, we 
suggest the following BibTex:

::

 @misc{GDT-Core,
       author = {Adam Goldstein and William H. Cleveland and Daniel Kocevski},
       title = {Gamma-ray Data Tools Core Package: v2.0.0},
       year = 2023,
       url = {https://github.com/USRA-STI/gdt-core}
 }

.. rubric:: Acknowledgments

The Gamma-ray Data Tools are partially funded through the NASA ADAP Grant 
80NSSC21K0651 and the NASA SMD Open Source Tools, Frameworks, and Libraries 
Grant 80NSSC22K1741.

***************
Getting Started
***************
.. toctree::
   :maxdepth: 1

   install

******************
User Documentation
******************

Data Types and Utilities
------------------------
.. toctree::
   :maxdepth: 1

   core/primitives
   core/collection
   core/healpix
   core/pha
   core/phaii
   core/response
   core/tte

File Access and Definitions
---------------------------
.. toctree::
   :maxdepth: 1

   core/file
   core/headers

Mission and Instrument Definitions
----------------------------------
.. toctree::
   :maxdepth: 1

   core/coords
   core/detectors
   core/geomagnetic
   core/heasarc
   core/time

Reduction and Analysis
----------------------
.. toctree::
   :maxdepth: 1
   
   core/background/background
   core/binning/binning   
   core/simulate/simulate
   core/spectra/spectra
   core/temporal/temporal
   
Plotting
--------
.. toctree::
   :maxdepth: 2
   
   core/plot/plot

*******
License
*******
.. toctree::
   :maxdepth: 1
   
   license


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
