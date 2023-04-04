.. _install:


Installation
============

..  Note:: Requires: Python >=3.8
            
           Tested on:
           
           * macOS Monterey (12.6.4) - Ventura (13.2.1)

How to Install
--------------

The GDT Core package can be installed from its tarfile gdt-core.tar.gz.

To install::

    $ pip install gdt-core.tar.gz

or to include requirements to build documentation::

    $ pip install gdt-core.tar.gz[docs]


Development
-----------
For development, we recommend the following, including installing the extra
dependencies for compiling the documentation and running unit tests::

    $ tar -xvzf gdt-core.tar.gz
    $ cd gdt-core
    $ pip install -e ".[docs,test]"

If you want to test your development, note that there are a number of data files
that are used in testing.  See :ref:`Downloading Test/Tutorial Data<download_test_data>` 
for downloading the test data.

To test your development, navigate to the gdt-core root directory and execute ::

    $ pytest


Documentation 
-------------
If installed using the ``docs`` option, the documentation can be compiled by
navigating to the gdt-core/docs directory and

    $ make html


.. _download_test_data:

Downloading Test/Tutorial Data
------------------------------
To download the data files used in the documentation and for testing, you need
to run the ``gdt-download-data`` script after installation. The downloader
script is designed so that you can download data from specific missions, or 
download all of the test/tutorial data.  To see the list of available missions::

    $ gdt-download-data --help

If you want to download the Fermi GBM test/tutorial data only, for example::

    $ gdt-download-data -m fermi-gbm

Or to download all of the data::
    
    $ gdt-download-data --all

The data are downloaded to a default directory. To access the data from the GDT, 
there is a variable at the main level that stores the path dictionary for each 
mission.  To access the Fermi GBM test data directory:

    >>> from gdt import test_data
    >>> gbm_path = test_data['fermi-gbm']
    
----

Quickstart
----------
To load the GDT Core package within your python environment, simply::
    
    import gdt.core


How to Uninstall
----------------

To uninstall::

    $ pip3 uninstall gdt-core

There are also a number of files for the tools that are copied into your 
``$HOME/.gammaray_data_tools`` directory.  You can delete these files if you 
wish.


Known Issues
------------
* **There appears to be some differences arising between installations on Mac ARM 
  processors (M1 and M2 chips) and other Mac or Linux processors.** As of now, 
  this only shows up when using some of the minimizers provided through 
  scipy.optimize.minimize for spectral fitting. Users can test for the presence
  of these differences by running the unit tests.  The known failures on Mac ARM
  processors are:
  
  * test_fitting.py::TestSpectralFitterOne::test_hessian
  * test_fitting.py::TestSpectralFitterOne::test_jacobian
  * test_fitting.py::TestSpectralFitterOne::test_residuals
  
  The current understanding is that differences arise in spectral fit values
  above machine precision, but represent < 1% relative errors on the fit values
  themselves. The exact origin of these differences is unclear, but may be 
  related to the underlying C or FORTRAN libraries and compilers that are used
  to compile scipy. Further investigation is ongoing.

* **The virtual environment is using your system ipython (or other package) 
  install.**  This can sometimes happen if you didn't install ipython (or other
  package) in the virtual environment.  Try installing ipython (or other package) 
  and restart your virtual environment.

* **You observe the following error**::
    
    ImportError: No module named '_tkinter'
  
  This is a situation where Matplotlib is using the ``tkinter`` backend for
  plotting.  You would see this error if you don't have ``tkinter`` installed. 
  You don't need to install ``tkinter`` if you don't want to; instead, you can
  create a file named `matplotlibrc` in your working directory that contains the
  following::
    
    backend : Agg