.. _spectra-parameters:
.. |Parameter| replace:: :class:`~gdt.core.data_primitives.Parameter`
.. |PhotonFlux| replace:: :class:`~gdt.core.spectra.parameters.PhotonFlux`
.. |PhotonFluence| replace:: :class:`~gdt.core.spectra.parameters.PhotonFluence`
.. |EnergyFlux| replace:: :class:`~gdt.core.spectra.parameters.EnergyFlux`
.. |EnergyFluence| replace:: :class:`~gdt.core.spectra.parameters.EnergyFluence`
.. |DetectorData| replace:: :class:`~gdt.core.spectra.parameters.DetectorData`
.. |ModelFit| replace:: :class:`~gdt.core.spectra.parameters.ModelFit`


*********************************************************
Spectral Parameters (:mod:`~gdt.core.spectra.parameters`)
*********************************************************
This module contains classes that support recording and containing spectral
parameter information resulting from spectral fits.

Fluxes and Fluences
===================
There are four specialized classes that inherit from |Parameter|: |PhotonFlux|,
|PhotonFluence|, |EnergyFlux|, and |EnergyFluence|. These classes provide some
small special functionality compared to the Parameter base class.  For example,
when we create a PhotonFlux, we must specify the energy range over which the
PhotonFlux is valid:

    >>> from gdt.core.spectra.parameters import PhotonFlux
    >>> # flux of 5 +/- 1 photons/cm^2/s over 50-300 keV
    >>> pflux = PhotonFlux(5.0, 1.0, (50.0, 300.0))
    >>> print(pflux)
    Photon Flux: 5.00 +/- 1.00 ph/cm^2/s

Note that the name, units, and support are already defined:

    >>> pflux.units
    'ph/cm^2/s'
    >>> pflux.support
    (0.0, inf)
    >>> pflux.energy_range
    (50.0, 300.0)
    
The same applies to the other flux and fluence parameters.

Detector Fit Info
=================
For a spectral fit, it is generally useful to collect metadata about the 
characteristics fo the detector and data used in the fit.  To that end,
|DetectorData| is primarily a container class that bundles together this
information.  To create a DetectorData object, we start out with the instrument
name, detector name, data type, and number of energy channels.

    >>> from gdt.core.spectra.parameters import DetectorData
    >>> det_data = DetectorData.from_data('Fermi, GBM', 'NAI_00', 'TTE', 128)
    >>> det_data
    <DetectorData: NAI_00; TTE;
     time range: (None, None) s;
     energy range: (None, None) keV>
    
While all info can be added to the DetectorData object on creation, the info
can also be added after initialization.  For example, we can update the time
and energy range:

    >>> det_data.time_range = (0.0, 10.0)
    >>> det_data.energy_range = (8.0, 900.0)
    >>> det_data
    
Several other attributes can be set such as the channel mask (which channels
were used for fitting), energy edges, response file that was used, and the best 
fit photon model counts and associated uncertainties.  Although not included in
the base DetectorData class, this class is useful for collecting this 
information so that it can be written to a (e.g. FITS) file.


Model Fit Info
==============
Similar to DetectorData, the |ModelFit| class serves as a container for the
info associated with a single spectral fit. It bundles together the fitted 
parameters and associated uncertainties, resulting fluxes and fluences, the
covariance matrix, and fit statistic and degrees-of-freedom.  To create a 
ModelFit object, we have to initialize with at least the name of the model and
the time range over which the data was fit:

    >>> from gdt.core.spectra.parameters import ModelFit
    >>> # fit a power law function from 0-10 s
    >>> model_fit = ModelFit.from_data('PowerLaw', (0.0, 10.0))
    >>> model_fit
    <ModelFit: PowerLaw>

Now we can add other pieces of information to the object, such as the parameter
fit values:

    >>> from gdt.core.data_primitives import Parameter
    >>> amp = Parameter(0.01, 2e-3, name='A', units='ph/cm^2/s/keV')
    >>> index = Parameter(-1.7, 0.3, name='index')
    >>> model_fit.parameters = [amp, index]
    >>> print(model_fit)
    PowerLaw
       A: 1.00e-02 +/- 2.00e-03 ph/cm^2/s/keV
       index: -1.70e+00 +/- 3.00e-01
       
    >>> model_fit.parameters
    [<Parameter: A>, <Parameter: index>]
       
We can also add flux information:

    >>> model_fit.photon_flux = PhotonFlux(1.659, 0.071, (50.0, 300.0))
    
And the fit statistic, etc:

    >>> model_fit.stat_name = 'PG-stat'
    >>> model_fit.stat_value = 143.45
    >>> model_fit.dof = 121

Reference/API
=============

.. automodapi:: gdt.core.spectra.parameters
   :inherited-members:


