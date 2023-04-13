======================================
Gamma-ray Data Tools - Core Components
======================================

The Gamma-ray Data Tools (GDT) is centralized toolkit for hard X-ray and
gamma-ray astrophysics data analysis, with a focus on providing a uniform
interface to the data provided by several different missions and instruments.

The GDT Core Package (``astro-gdt``) contains the core components of the GDT that
can be utilized for various instruments and is a generalized version of the
`Fermi GBM Data Tools <https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs>`_.
Individual mission or instrument packages can be developed using ``astro-gdt``
and released under the ``gdt.missions`` namespace (see ``astro-gdt-fermi`` as an example).


Normal Installation
-------------------

If you don't plan to contribute code to the project, the recommended install method is installing from PyPI using:

.. code-block:: sh

   pip install astro-gdt
   gdt-data init

The ``gdt-data init`` is required to initialize the library after installation.


Contributing Code or Documentation
----------------------------------

If you plan to help with the development or documentation of astro-gdt, then please visit our github site at
https://github.com/USRA-STI/gdt-core.
