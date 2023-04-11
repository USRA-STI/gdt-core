.. _core-file:
.. |FitsFileContextManager| replace:: :class:`~gdt.core.file.FitsFileContextManager`
.. |HDUList| replace:: HDUList
.. _HDUList: https://docs.astropy.org/en/stable/io/fits/api/hdulists.html

**************************************
The File Module (:mod:`gdt.core.file`)
**************************************

Introduction
============
The File module contains, of primary importance, the |FitsFileContextManager|
abstract class. As the name suggests, it is a context manager for accessing 
FITS files. While this class is intended to be inherited, we can demonstrate
how this class works and is intended to be used.

Examples
--------
Let's assume we have the FITS file "glg_ctime_nb_bn120415958_v00.pha" in our
current directory.  You can access it with the following:

    >>> from gdt.core import data_path
    >>> from gdt.core.file import FitsFileContextManager
    >>> filepath = data_path.joinpath('fermi-gbm').joinpath("glg_ctime_nb_bn120415958_v00.pha")
    >>> f = FitsFileContextManager.open(filepath)
    >>> f
    <FitsFileContextManager(filename="glg_ctime_nb_bn120415958_v00.pha") at 0x109d25410>

The object contains some convenience functions:

    >>> # number of HDUs
    >>> f.num_hdus
    4
    
    >>> # get the column names in the second HDU (zero-indexed)
    >>> f.get_column_names(2)
    ('COUNTS', 'EXPOSURE', 'QUALITY', 'TIME', 'ENDTIME')
    
    >>> # retrieve the 'EXPOSURE' column
    >>> f.column(2, 'EXPOSURE')
    array([0.25506625, 0.25459924, 0.25496757, ..., 0.2551209 , 0.25507143,
           0.25600004], dtype=float32)
    
    >>> # retrieve both the 'TIME' and 'ENDTIME' columns in a single array
    >>> f.columns_as_array(2, ['TIME', 'ENDTIME'])
    array([[3.56222662e+08, 3.56222662e+08],
           [3.56222662e+08, 3.56222662e+08],
           [3.56222662e+08, 3.56222663e+08],
           ...,
           [3.56224561e+08, 3.56224561e+08],
           [3.56224561e+08, 3.56224562e+08],
           [3.56224562e+08, 3.56224562e+08]])

The underlying Astropy |HDUList|_ can be accessed directly:

    >>> f.hdulist
    <astropy.io.fits.hdu.image.PrimaryHDU object at 0x103d31dd0>, 
    <astropy.io.fits.hdu.table.BinTableHDU object at 0x108b294d0>, 
    <astropy.io.fits.hdu.table.BinTableHDU object at 0x108af4810>, 
    <astropy.io.fits.hdu.table.BinTableHDU object at 0x108b300d0>]
    
Whenever we access the file this way, we must remember to close it:

    >>> f.close()

The power behind using a context manager is that we can open the file in a 
``with`` block and the file is closed properly when we exit the block or hit an
error:

    >>> with FitsFileContextManager(filepath) as f:
    >>>     times = f.column(2, 'TIME')
    >>>
    >>> print(times)
    array([3.56222662e+08, 3.56222662e+08, 3.56222662e+08, ...,
           3.56224561e+08, 3.56224561e+08, 3.56224562e+08])


Reference/API
=============

.. automodapi:: gdt.core.file
   :inherited-members:

