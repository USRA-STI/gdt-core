.. _install:

Installation
============

..  Note:: Requires: Python >=3.7
            
           Tested on:
           
           * macOS High Sierra (10.12.0) - Ventura (13.2.1)
           
           * Ubuntu 16.04 - 18.04
           
           * Windows 10 Subsystem for Linux


How to Install
--------------

The GDT can be installed from its tarfile gammaray_data_tools-2.0.0.tar.gz.

To install::

    $ pip install gammaray_data_tools-2.0.0.tar.gz

or to include requirments to build documentation::

    $ pip install gammaray_data_tools-2.0.0.tar.gz[docs]

For development, we recommend the following::

    $ tar -xvzf gammaray_data_tools-2.0.0.tar.gz
    $ cd gdt-2.0.0
    $ pip install -e ".[docs]"


If you want to test your development, note that there are a number of data files
that are used in testing.  See :ref:`Downloading Test/Tutorial Data<download_test_data>` 
for downloading the test data.

To test your development::

    $ pytest test

If you installed the documentation requirements (with [docs]), then you can go 
ahead and install the documentation and tutorial within your home directory with::

    $ install-gdt-docs

----

How to Uninstall
----------------

To uninstall::

    $ pip3 uninstall gdt

There are also a number of files (documentation, notebooks, etc.) for the tools
that are copied into your ``$HOME/.gammaray_data_tools`` directory.  You can 
delete these files if you wish.

Documentation 
-------------
On successful installation, you can launch the local HTML documentation by
calling::

    $ gdt-docs

from the command line.

Some of the documentation uses real data files to demonstrate functionality. 
See the next section on downloading this data.


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
To load the GDT within your python environment, simply::
    
    import gdt
    

Known Issues
------------
* **When running a notebook in Linux, you observe a similar error**::
    
    %matplotlib notebook                                                               
    Warning: Cannot change to a different GUI toolkit: notebook. Using osx 
    instead. This is due to some backend plotting issue with Jupyter notebook 
    on Linux. Remove the ``%matplotlib inline`` in the notebook cell and 
    re-evaluate the cell *twice* to see the plot.


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