# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
# Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland and Adam Goldstein
#               Universities Space Research Association
#               Science and Technology Institute
#               https://sti.usra.edu
#
# Developed by: Daniel Kocevski
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#
import numpy as np
from gdt.core.data_primitives import Parameter
from gdt.core.types import set_by_dict

__all__ = ['DetectorData', 'EnergyFlux', 'EnergyFluence', 'ModelFit',
           'PhotonFluence', 'PhotonFlux']


class PhotonFlux(Parameter):
    """A photon flux class.
    
    Parameters:
        value (float): The central flux value
        uncert (float or 2-tuple): The 1-sigma uncertainty
        energy_range (tuple): A 2-tuple (low, high) for the energy range
    """

    def __init__(self, value, uncert, energy_range):
        super().__init__(value, uncert, name='Photon Flux', units='ph/cm^2/s',
                         support=(0.0, np.inf))
        self._energy_range = energy_range

    @property
    def energy_range(self):
        """(tuple): The enery range (low, high)"""
        return self._energy_range


class PhotonFluence(Parameter):
    """A photon fluence class.
    
    Parameters:
        value (float): The central fluence value
        uncert (float or 2-tuple): The 1-sigma uncertainty
        energy_range (tuple): A 2-tuple (low, high) for the energy range
    """

    def __init__(self, value, uncert, energy_range):
        super().__init__(value, uncert, name='Photon Fluence', units='ph/cm^2',
                         support=(0.0, np.inf))
        self._energy_range = energy_range

    @property
    def energy_range(self):
        """(tuple): A 2-tuple (low, high) for the energy range"""
        return self._energy_range


class EnergyFlux(Parameter):
    """An energy flux class.
    
    Parameters:
        value (float): The central flux value
        uncert (float or 2-tuple): The 1-sigma uncertainty
        energy_range (tuple): A 2-tuple (low, high) for the energy range
    """

    def __init__(self, value, uncert, energy_range):
        super().__init__(value, uncert, name='Energy Flux', units='erg/cm^2/s',
                         support=(0.0, np.inf))
        self._energy_range = energy_range

    @property
    def energy_range(self):
        """(tuple): A 2-tuple (low, high) for the energy range"""
        return self._energy_range


class EnergyFluence(Parameter):
    """An energy fluence class.
    
    Parameters:
        value (float): The central fluence value
        uncert (float or 2-tuple): The 1-sigma uncertainty
        energy_range (tuple): A 2-tuple (low, high) for the energy range
    """

    def __init__(self, value, uncert, energy_range):
        super().__init__(value, uncert, name='Energy Fluence', units='erg/cm^2',
                         support=(0.0, np.inf))
        self._energy_range = energy_range

    @property
    def energy_range(self):
        """(tuple): A 2-tuple (low, high) for the energy range"""
        return self._energy_range


class ModelFit:
    """A container for the info resulting from a spectral fit.
    """

    def __init__(self):
        self._name = None
        self._time_range = None

        self._parameters = []
        self._photon_flux = None
        self._energy_flux = None
        self._photon_fluence = None
        self._energy_fluence = None
        self._flux_energy_range = None
        self._stat_name = None
        self._stat_value = None
        self._dof = None
        self._covariance = None

    @property
    def covariance(self):
        """(np.array): The covariance matrix of the fit"""
        return self._covariance

    @covariance.setter
    def covariance(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('covariance must an array')
        if len(val.shape) != 2:
            raise ValueError('covariance must be a n x n array')
        if val.shape[0] != val.shape[1]:
            raise ValueError('covariance must be a n x n array')

        self._covariance = val

    @property
    def dof(self):
        """(int): The degrees-of-freedom of the fit"""
        return self._dof

    @dof.setter
    def dof(self, val):
        try:
            int_val = int(val)
        except(ValueError, TypeError):
            raise TypeError('dof must be an integer')
        self._dof = int_val

    @property
    def energy_fluence(self):
        """(:class:`EnergyFluence`): The energy fluence"""
        return self._energy_fluence

    @energy_fluence.setter
    def energy_fluence(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('energy_fluence must be of Parameter type')
        self._energy_fluence = val

    @property
    def energy_flux(self):
        """(:class:`EnergyFlux`): The energy flux"""
        return self._energy_flux

    @energy_flux.setter
    def energy_flux(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('energy_flux must be of Parameter type')
        self._energy_flux = val

    @property
    def flux_energy_range(self):
        """(tuple): The energy range of the flux and fluence, (low, high)"""
        return self._flux_energy_range

    @flux_energy_range.setter
    def flux_energy_range(self, val):
        if not isinstance(val, (list, tuple)):
            raise TypeError('flux_energy_range must be a 2-tuple')
        else:
            if len(val) != 2:
                raise ValueError('flux_energy_range must be a 2-tuple')
        self._flux_energy_range = val

    @property
    def name(self):
        """(str): The name of the model"""
        return self._name

    @property
    def parameters(self):
        """(list of :class:`~gdt.core.data_primitives.Parameter`): A list of 
        model parameters"""
        return self._parameters

    @parameters.setter
    def parameters(self, val):
        if not isinstance(val, (list, tuple)):
            raise TypeError('parameters must be a list of parameters')
        for p in val:
            if not isinstance(p, Parameter):
                raise TypeError('parameters must be of Parameter type')
        self._parameters = val

    @property
    def photon_fluence(self):
        """(:class:`PhotonFluence`): The photon fluence"""
        return self._photon_fluence

    @photon_fluence.setter
    def photon_fluence(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('photon_fluence must be of Parameter type')
        self._photon_fluence = val

    @property
    def photon_flux(self):
        """(:class:`PhotonFlux`): The photon flux"""
        return self._photon_flux

    @photon_flux.setter
    def photon_flux(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('photon_flux must be of Parameter type')
        self._photon_flux = val

    @property
    def stat_name(self):
        """(str): The name of the fit statistic"""
        return self._stat_name

    @stat_name.setter
    def stat_name(self, val):
        self._stat_name = str(val)

    @property
    def stat_value(self):
        """(float): The fit statistic value"""
        return self._stat_value

    @stat_value.setter
    def stat_value(self, val):
        try:
            float_val = float(val)
        except (ValueError, TypeError):
            raise TypeError('stat_value must be a float')
        self._stat_value = float_val

    @property
    def time_range(self):
        """(float, float): The time range of the model fit, (low, high)"""
        return self._time_range

    @classmethod
    def from_data(cls, name, time_range, **kwargs):
        """Create a ModelFit object from data.
        
        Args:
            name (str): The name of the model
            time_range (float, float): The time range of the model fit, 
                                       (low, high)
            parameters (list, optional): A list of model parameters
            photon_flux (:class:`PhotonFlux`, optional): The photon flux
            energy_flux (:class:`EnergyFlux`, optional): The energy flux
            photon_fluence (:class:`PhotonFluence`, optional): The photon fluence
            energy_fluence (:class:`EnergyFluence`, optional): The energy fluence
            flux_energy_range (tuple, optional): The energy range of the flux
                                                 and fluence, (low, high)
            stat_name (str, optional): The name of the fit statistic
            stat_value (float, optional): The fit statistic value
            dof (int, optional): The degrees-of-freedom of the fit
            covariance (np.array, optional): The covariance matrix of the fit
        
        Returns:
            (:class:`ModelFit`)
        """
        obj = cls()
        obj._name = str(name)
        if not isinstance(time_range, (list, tuple)):
            raise TypeError('time_range must be a 2-tuple')
        else:
            if len(time_range) != 2:
                raise ValueError('time_range must be a 2-tuple')
        obj._time_range = time_range

        set_by_dict(obj, kwargs, no_private=False)
        return obj

    def parameter_list(self):
        """Return the list of parameter names
        
        Returns:
            (list): The parameter names
        """
        return [param.name for param in self.parameters]

    def __repr__(self):
        s = '<{0}: {1}>'.format(self.__class__.__name__, self.name)
        return s

    def __str__(self):
        param_str = '\n   '.join([str(param) for param in self.parameters])
        return '{0}\n   {1}'.format(self.name, param_str)


class DetectorData:
    """A container for detector info used in a spectral fit.    
    """

    def __init__(self):
        self._instrument = None
        self._detector = None
        self._datatype = None
        self._filename = None
        self._numchans = None

        self._active = True
        self._response = ''
        self._time_range = (None, None)
        self._energy_range = (None, None)
        self._channel_range = (None, None)
        self._channel_mask = None
        self._energy_edges = None
        self._photon_counts = None
        self._photon_model = None
        self._photon_errors = None

    @property
    def active(self):
        """(bool, optional): True if the detector is used in the fit"""
        return self._active

    @active.setter
    def active(self, val):
        self._active = bool(val)

    @property
    def channel_mask(self):
        """(np.array): A Boolean array marking the channels that were used"""
        return self._channel_mask

    @channel_mask.setter
    def channel_mask(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('channel_mask must be an array')
        self._channel_mask = val.astype(bool)

    @property
    def channel_range(self):
        """(int, int): The energy channel range of the data"""
        return self._channel_range

    @channel_range.setter
    def channel_range(self, val):
        if not isinstance(val, (tuple, list)):
            raise TypeError('channel_range must be a 2-tuple')
        elif len(val) != 2:
            raise ValueError('channel_range must be a 2-tuple')
        elif val[0] > val[1]:
            raise ValueError('channel_range must be of form (low, high)')
        else:
            pass

        self._channel_range = val

    @property
    def datatype(self):
        """(str): The name of the datatype"""
        return self._datatype

    @property
    def detector(self):
        """(str): The name of the detector"""
        return self._detector

    @property
    def energy_edges(self):
        """(np.array): The edges of the energy channels"""
        return self._energy_edges

    @energy_edges.setter
    def energy_edges(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('energy_edges must be an array')
        self._energy_edges = val

    @property
    def energy_range(self):
        """(float, float): The energy range of the data used"""
        return self._energy_range

    @energy_range.setter
    def energy_range(self, val):
        if not isinstance(val, (tuple, list)):
            raise TypeError('energy_range must be a 2-tuple')
        elif len(val) != 2:
            raise ValueError('energy_range must be a 2-tuple')
        elif val[0] > val[1]:
            raise ValueError('energy_range must be of form (low, high)')
        else:
            pass

        self._energy_range = val

    @property
    def filename(self):
        """(str): The filename of the data file"""
        return self._filename

    @filename.setter
    def filename(self, val):
        self._filename = val

    @property
    def instrument(self):
        """(str): The name of the instrument"""
        return self._instrument

    @property
    def num_chans(self):
        """(int): Number of energy channels used"""
        return self._numchans

    @property
    def photon_counts(self):
        """(np.array): The deconvolved photon counts for the detector"""
        return self._photon_counts

    @photon_counts.setter
    def photon_counts(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('photon_counts must be an array')
        self._photon_counts = val

    @property
    def photon_errors(self):
        """(np.array): The deconvolved photon count errors for the detector"""
        return self._photon_errors

    @photon_errors.setter
    def photon_errors(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('photon_errors must be an array')
        self._photon_errors = val

    @property
    def photon_model(self):
        """(np.array): The photon model for the detector"""
        return self._photon_model

    @photon_model.setter
    def photon_model(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('photon_model must be an array')
        self._photon_model = val

    @property
    def response(self):
        """(str): The filename of the detector response"""
        return self._response

    @response.setter
    def response(self, val):
        self._response = str(val)

    @property
    def time_range(self):
        """(float, float): The time range of the data used"""
        return self._time_range

    @time_range.setter
    def time_range(self, val):
        if not isinstance(val, (tuple, list)):
            raise TypeError('time_range must be a 2-tuple')
        elif len(val) != 2:
            raise ValueError('time_range must be a 2-tuple')
        elif val[0] > val[1]:
            raise ValueError('time_range must be of form (low, high)')
        else:
            pass

        self._time_range = val

    @classmethod
    def from_data(cls, instrument, detector, datatype, numchans, **kwargs):
        """Create a DetectorData object from data.
        
        Args:
            instrument (str): The name of the instrument
            detector (str): The name of the detector
            datatype (str): The name of the datatype
            filename (str): The filename of the data file
            numchans (int): Number of energy channels used
            active (bool, optional): True if the detector is used in the fit
            channel_mask (np.array, optional): A Boolean array marking the 
                                               channels that were used
            channel_range (tuple, optional): The energy channel range of the 
                                             data
            energy_edges (np.array, optional): The edges of the energy channels
            energy_range (tuple, optional): The energy range of the data used
            photon_counts (np.array, optional): The deconvolved photon counts 
                                                for the detector
            photon_errors (np.array, optional): The deconvolved photon count 
                                                errors for the detector
            photon_model (np.array, optional): The photon model for the detector
            response (str, optional): The filename of the detector response
            time_range (tuple, optional): The time range of the data used
        
        Returns:
            (:class:`DetectorData`)
        """
        obj = cls()
        obj._instrument = str(instrument)
        obj._detector = str(detector)
        obj._datatype = str(datatype)
        obj._numchans = int(numchans)
        set_by_dict(obj, kwargs, no_private=False)
        return obj

    def __repr__(self):
        s = '<{0}: {1}; {2};\n'.format(self.__class__.__name__, self.detector,
                                       self.datatype)
        s += ' time range: {} s;\n'.format(self.time_range)
        s += ' energy range: {} keV>'.format(self.energy_range)
        return s
