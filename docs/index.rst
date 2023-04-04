.. The Gamma-ray Data Tools documentation master file, created by
   sphinx-quickstart on Sun Mar 27 17:59:23 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to documentation for the Gamma-ray Data Tools Core Package!
===================================================================
The Gamma-ray Data Tools (GDT) is centralized toolkit for hard X-ray and 
gamma-ray astrophysics data analysis, with a focus on providing a uniform 
interface to the data provided by several different missions and instruments.

The GDT Core Package (``gdt-core``) contains the core components of the GDT that
can be utilized for various instruments.  Individual mission or instrument
packages can be developed using the ``gdt-core`` and released under the ``gdt``
namespace (see ``gdt-fermi`` as an example).

The documentation linked below walks through all of sub-packages and modules
within the gdt-core, and developers should take special note of the 
"For Developers" sections that detail how to subclass or design plugin functions
or classes for use in new instrument packages.


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
   
Plotting
--------
.. toctree::
   :maxdepth: 2
   
   core/plot/plot

******************
License
******************
.. toctree::
   :maxdepth: 1
   
   license


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
