======================================
Gamma-ray Data Tools - Core Components
======================================

This software is not subject to EAR.

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


Contributing to the GDT
-----------------------

Community contributions to the GDT are welcome and make the GDT a better
toolkit for everyone to use.  There are various types of contributions.  These
include:

* Bug fixes
* New features
* Documentation improvements
* API improvements
* New mission packages

For each contribution, we consider it best practice to first open an issue (with
an appropriate ``label``).  This may be as simple as notifying us (and other
users!) that there is a bug.  Perhaps you also have a proposed solution to
fixing the issue.  If so, then you can outline your proposed fix when you open
the issue, and we will give you feedback about whether we think that would be a
useful and appropriate fix.

**Simple typographical fixes and clarifications made in the documentation do not require
the creation of an issue.**

If your proposed modifications are significant (e.g. sizable document change, API improvements,
new feature), we highly recommend that you detail the propose change in the issue and wait for feedback.
This is to save you precious time in the event that we decide not to accept your proposed solution.
Often the proposed solution can break functionality elsewhere or can be simplified, and we would like to have
a chance to provide useful feedback before you begin coding. Waiting for feedback
is not a requirement, but merely reduces the chance of additional changes prior to having
your pull request accepted.

If you are submitting code modifications, we require that you create a unit test to confirm
expected operation of the code if those modifications aren't already covered by an
existing unit test.

The usual sequence of events are:

1. Create an issue describing the proposed changes.
2. Waiting for feedback if desired.
3. Create a fork from main branch.
4. Use your fork to add your changes to the code base.
5. Create unit tests showing that your changes work as intended.
6. Create `Pull Request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`_ with a comment explaining how it closes the issue you created.



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

