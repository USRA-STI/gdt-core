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

The full documentation can be found `here <https://astro-gdt.readthedocs.io/en/latest/>`_.

Normal Installation
-------------------

If you don't plan to contribute code to the project, the recommended install method is installing from PyPI using:

.. code-block:: sh

   pip install astro-gdt
   gdt-data init

The ``gdt-data init`` is required to initialize the library after installation.


Setting up a development environment
------------------------------------

If you do want to contribute code to this project (and astro-gdt), you can use the following commands to quickly setup a
development environment:

.. code-block:: sh

   mkdir gdt-devel
   cd gdt-devel
   python -m venv venv
   . venv/bin/activate
   pip install --upgrade pip setuptools wheel
   git clone git@github.com:USRA-STI/gdt-core.git
   pip install -e gdt-core/
   gdt-data init
   pip install -r gdt-core/requirements.txt

This should result in git-devel having the following directory structure::

   .
   ├── venv
   └── gdt-core

with gdt-core installed in the virtual environment named venv.

Writing Extensions using Namespace Packaging
--------------------------------------------
Gamma-ray Data Tools encourages missions to write extensions using namespace packages. Please use our
`Fermi extension <https://github.com/USRA-STI/gdt-fermi>`_ as an example of how we expect other missions to contribute
extensions to the Gamma-ray Data Tools.

The extension package should contain a directory 'gdt' with a subdirectory 'missions' which will hold the extension code
in a package directory named after the mission.

For example, GDT-Fermi has the following directory layout::

  .
  ├── config
  ├── dist
  ├── docs
  ├── src
  │   └── gdt
  │      └── missions
  │          └── fermi
  │              ├── gbm
  │              │   └── __init__.py
  │              ├── lat
  │              │   └── __init__.py
  │              └── __init__.py
  └── tests
    └── missions
        └── fermi


Since GDT-Fermi uses namespace packaging, both ``src/gdt`` and  ``src/gdt/missions`` do not contain a file named
``__init__.py``. This is because they are Namespace packages.

Notice that directory ``src/gdt/mission/fermi`` and its subdirectories contains an `__init__.py` file
signalling to Python that those directories are regular packages.

You can learn more about Namespace packages by reading `PEP-420 <https://peps.python.org/pep-0420/>`_.

Helping with Documentation
--------------------------

You can contribute additions and changes to the documentation. In order to use sphinx to compile the documentation
source files, we recommend that you install the packages contained within ``requirments.txt``.

To compile the documentation, use the following commands:

.. code-block:: sh

   cd gdt-core/docs
   make html

