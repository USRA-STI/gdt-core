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
from abc import ABC, abstractmethod
from typing import List, Tuple

import numpy as np
from scipy.integrate import trapz, quad
from scipy.special import erfc

__all__ = ['Band', 'BandOld', 'BlackBody', 'BrokenPowerLaw', 'Comptonized',
           'DoubleBrokenPowerLaw', 'GaussLine', 'GaussLineMult',
           'GaussianLog', 'GaussianLogVaryingFWHM', 'HighEnergyCutoff',
           'LogNormal', 'LorentzLineMult', 'LowEnergyCutoff', 'OTTB', 'OTTS',
           'PowerLaw', 'PowerLawMult', 'SmoothlyBrokenPowerLaw',
           'SunyaevTitarchuk', 'TanakaPulsar', 'YangSoongPulsar', 'Function',
           'SuperFunction']


class Function(ABC):
    """Abstract class for a 1-D function (e.g. photon model).
    
    Sub-classes should be designed with the following requirements:
    
    #. Attribute :attr:`nparams` containing the total number of parameters of the
       function (even parameters that are fixed to a particular value during fitting)
    #. Method :meth:`eval` that accepts two arguments: an array of parameter values
       of length :attr:`nparams`; and the array of abscissa values.
       
    Subclasses can optionally add the following for complete and enhanced
    functionality:
    
    #. Attribute :attr:`param_list` that is a list of parameter names in the order
       they are accepted in the :meth:`eval` method. Each entry in :attr:`param_list` is
       a tuple of the form (<param_name>, <units>, <description>).  If no units
       or description are needed, those fields should be a null string ''. If
       omitted, these will be filled to default values ('param0', '', ''), etc.
    #. Attribute :attr:`default_values` that is a list of default parameter values. 
       This is primarily used for parameter estimation. Default is 1.0 for each
       parameter
    #. Attribute :attr:`delta_abs` that is a list of max absolute fitting step for 
       each parameter.  This is used for parameter estimation. Default is 0.1
       for each parameter.
    #. Attribute :attr:`delta_rel` that is a list of max relative fitting step for
       each parameter.  This is used for paramter estimation.  Default is 0.01
       for each parameter.
    #. Attribute :attr:`min_values` that is a list of minimum values each parameter
       can take.  This is used primarily for parameter estimation. Default is
       -np.inf for each parameter.
    #. Attribute :attr:`max_values` that is a list of maximum values each parameter
       can take.  This is used primarily for parameter estimation. Default is
       np.inf for each parameter.
    #. Attribute :attr:`free` that is a boolean list indicating if each parameter
       is allowed to be free for fitting.  False indicates that the parameter
       will be fixed at its default value. This is used for parameter estimation. 
       Default is True for each parameter.
    """

    #: Number of parameters in the function.
    nparams: int = None

    #: The list of parameters.  Each element is of form (parameter name, units, short description)
    param_list: List[Tuple[str, str, str]] = None

    #: The name of the function.
    name: str = None

    #: The default parameter values. If not defined, each defaults to 1.0
    default_values: List[float] = None

    #: The maximum absolute fitting step size for each parameter. If not defined, each default to 0.1
    delta_abs: List[float] = None

    #: The maximum relative fitting step size for each parameter. If not defined, each default to 0.01
    delta_rel: List[float] = None

    #: The minimum possible value for each parameter. If not defined, each default to -np.inf
    min_values: List[float] = None

    #: The maximum possible value for each parameter. If not defined, each default to np.inf
    max_values: List[float] = None

    #: For each parameter, mark True if left free for fitting, otherwise mark False to fix at the default
    #: value. If not defined, each default to True."""
    free: List[bool] = None

    def __init__(self):
        self._num_components = 1

        if self.nparams is None:
            raise AttributeError('Function MUST have nparams defined')

        if self.param_list is None:
            self.param_list = [('param' + str(i), '', '') for i in range(self.nparams)]

        if self.name is None:
            self.name = self.__class__.__name__

        if self.default_values is None:
            self.default_values = [1.0] * self.nparams

        if self.delta_abs is None:
            self.delta_abs = [0.1] * self.nparams

        if self.delta_rel is None:
            self.delta_rel = [0.01] * self.nparams

        if self.min_values is None:
            self.min_values = [-np.inf] * self.nparams

        if self.max_values is None:
            self.max_values = [np.inf] * self.nparams

        if self.free is None:
            self.free = [True] * self.nparams

    @property
    def num_components(self):
        """(int): Number of function components. For a normal function this is 
                  one, for a :class:`SuperFunction`, this can be > 1."""
        return self._num_components

    @abstractmethod
    def eval(self, params, x):  # pragma: no cover
        """Evaluate the function.
        
        Args:
            params (iterable): The parameter vector of the free parameters 
            x (np.array): The values at which to evaluate the function 
        
        Returns:        
            (np.array)
        """
        pass

    def fit_eval(self, params, x, **kwargs):
        """Evaluate the function for a fit.  This differs from :meth:`eval` 
        because the ``params`` vector should only contain the free fit 
        parameters. The parameters that are set to fixed will automatically be 
        passed to :meth:`eval` at their fixed value.  This enables a fit to be 
        performed with fixed parameters and the associated jacobian and 
        covariance matrices will have the appropriate shape.
        
        Args:
            params (iterable):  The parameter vector of the free parameters 
            x (np.array): The values at which to evaluate the function 
        
        Returns:        
            (np.array)
        """
        full_params = np.zeros(self.nparams)
        free = np.array(self.free)
        full_params[free] = params
        full_params[~free] = np.array(self.default_values)[~free]
        return self.eval(full_params, x, **kwargs)

    def integrate(self, params, erange, num_points=1000, log=True,
                  energy=False):
        """Integrate the function over energy to produce a flux.
        
        Args:
            params (iterable): The parameter vector, which can be of length 
                               :attr:`nparams` or equal to the length of the 
                               non-fixed (free) parameters.
            erange (float, float): The energy range over which to integrate
            num_points (int, optional): The number of points used for 
                                        integration. Default is 1000.
            log (bool, optional): If True, the integration grid is made in 
                                  log-space, otherwise it is linear. Default 
                                  is True.
            energy (bool, optional): If True, perform integration of E*F(E) to 
                                     produce an energy flux, otherwise 
                                     integration is over F(E) to produce a 
                                     photon flux. Default is False.
        
        Returns:        
            (float)
        """
        # If given all params, use eval
        if len(params) == self.nparams:
            func = self.eval
        # otherwise, use fit_eval (auto apply fixed params)
        elif len(params) == np.sum(self.free):
            func = self.fit_eval
        else:
            raise ValueError('Incorrect number of parameters')

        # evenly-spaced energies in log space
        if log:
            energies = np.logspace(*np.log10(erange), num_points)
        # evenly-spaced energies in linear space
        else:
            energies = np.linspace(*erange, num_points)

        spectrum = func(params, energies)

        # photon flux
        if not energy:
            flux = trapz(spectrum, energies)
        # energy flux
        else:
            flux = trapz(energies * spectrum, energies) * 1.602e-9

        return flux

    def parameter_bounds(self, apply_state=True):
        """The valid parameter bounds in the form [(min, max), ...]
        
        Args:
            apply_state (bool, optional): 
                If True, will only return the bounds for free parameters. 
                Default is True.
        
        Returns:
            ([(float, float), ...]):
        """
        bounds = list(zip(*(self.min_values, self.max_values)))
        if apply_state:
            bounds = [bounds[i] for i in range(self.nparams) if self.free[i]]
        return bounds

    @classmethod
    def add(cls, func1, func2):
        """Add two Functions. Can also use the ``+`` operator.
        
        Args:
            func1 (:class:`Function`): A function 
            func2 (:class:`Function`): A function to add to ``func1``.
        
        Returns:
            (:class:`SuperFunction`)
        """
        obj = cls._super_function(func1, func2, np.sum)
        obj.names = [func1.name, func2.name]
        obj.name = ' + '.join(obj.names)
        return obj

    @classmethod
    def multiply(cls, func1, func2):
        """Multiply two Functions. Can also use the ``*`` operator.
        
        Args:
            func1 (:class:`Function`): A function 
            func2 (:class:`Function`): A function to multiply to ``func1``.
        
        Returns:
            (:class:`SuperFunction`)
        """
        obj = cls._super_function(func1, func2, np.prod)
        obj.names = [func1.name, func2.name]
        obj.name = ' * '.join(obj.names)
        return obj

    @classmethod
    def _super_function(cls, func1, func2, operator):
        """Creates a SuperFunction from two functions (one may itself be
        a SuperFunction) and an operator.
        
        Args:
            func1 (:class:`Function` or :class:`SuperFunction`): A function 
            func2 (:class:`Function` or :class:`SuperFunction`): Another function.
            operator (<function>): The operator to be performed on the two functions.
        
        Returns:
            :class:`SuperFunction`: A function containing the product of the \
                                    two functions
        """
        obj = SuperFunction()
        obj._num_components = func1._num_components + func2._num_components
        obj.nparams = func1.nparams + func2.nparams

        # func1 could be a SuperFunction or just a simple Function        
        if isinstance(func1, SuperFunction):
            obj._operator = func1._operator
            obj._sub_names = func1._sub_names
            obj._sub_param_list = func1._sub_param_list
            obj._sub_nparams = func1._sub_nparams
            obj._sub_functions = func1._sub_functions
            obj._operator.append(operator)
            obj._sub_names.append(func2.name)
            obj._sub_param_list.append(func2.param_list)
            obj._sub_nparams.append(func2.nparams)
            obj._sub_functions.append(func2.eval)
        elif isinstance(func2, SuperFunction):
            obj._operator = [operator]
            obj._sub_names = [func1.name]
            obj._sub_param_list = [func1.param_list]
            obj._sub_nparams = [func1.nparams]
            obj._sub_functions = [func1.eval]
            obj._operator.extend(func2._operator)
            obj._sub_names.extend(func2._sub_names)
            obj._sub_param_list.extend(func2._sub_param_list)
            obj._sub_nparams.extend(func2._sub_nparams)
            obj._sub_functions.extend(func2._sub_functions)
        else:
            obj._operator = [operator]
            obj._sub_names = [func1.name]
            obj._sub_param_list = [func1.param_list]
            obj._sub_nparams = [func1.nparams]
            obj._sub_functions = [func1.eval]
            obj._sub_names.append(func2.name)
            obj._sub_param_list.append(func2.param_list)
            obj._sub_nparams.append(func2.nparams)
            obj._sub_functions.append(func2.eval)

        # modify param lists to have component names
        obj.param_list = obj._modify_param_list()

        obj.default_values = func1.default_values + func2.default_values
        obj.delta_abs = func1.delta_abs + func2.delta_abs
        obj.delta_rel = func1.delta_rel + func2.delta_rel
        obj.min_values = func1.min_values + func2.min_values
        obj.max_values = func1.max_values + func2.max_values
        obj.free = func1.free + func2.free
        return obj

    def __add__(self, other):
        return Function.add(self, other)

    def __mul__(self, other):
        return Function.multiply(self, other)

    def __repr__(self):
        s = "<{0}: {1} parameters;\n".format(self.name, self.nparams)
        s += " Defaults: "
        for i in range(self.nparams):
            if i == 0:
                spaces = ''
            else:
                spaces = ' ' * 11
            if i == self.nparams - 1:
                end = '>'
            else:
                end = ';\n'
            start = '*' if not self.free[i] else ''
            s += "{0}{1}{2} = {3} {4}{5}".format(spaces, start,
                                                 self.param_list[i][0],
                                                 self.default_values[i],
                                                 self.param_list[i][1], end)
        return s


class SuperFunction(Function):
    """Abstract class for a sum or product of functions.
    """
    nparams = 0

    def __init__(self):
        super().__init__()
        self._sub_names = []
        self._sub_nparams = []
        self._sub_param_list = []
        self._sub_functions = []
        self._operator = []

    def eval(self, params, x, components=False):
        """Evaluate the function. Can return the evaluation of the overall 
        function or the evaluation of each component to the function.
        
        Args:
            params (iterable): The parameter vector of the free parameters 
            x (np.array): The values at which to evaluate the function 
            components (bool, optional): If True, evaluate for each component. 
                                         Default is False.
        
        Returns:        
            (np.array or list of np.array)
        """
        ind = np.concatenate(([0], np.cumsum(self._sub_nparams)))
        evals = [self._sub_functions[i](params[ind[i]:ind[i + 1]], x) for i in range(len(self._sub_functions))]
        if components:
            return evals
        else:
            the_eval = evals[0]
            for i in range(self.num_components - 1):
                the_eval = self._operator[i]((the_eval, evals[i + 1]), axis=0)
            return the_eval

    def _modify_param_list(self):
        """Modify the parameter list so the parmeter names contain the parent
        function name.
        
        Returns: 
            [(str, str, str), ...]: The modified parameter list
        """
        param_list = []
        for name, plist in zip(self._sub_names, self._sub_param_list):
            new_list = [(name + ': ' + p[0], p[1], p[2]) for p in plist]
            param_list.extend(new_list)
        return param_list


# ---- functions
class PowerLaw(Function):
    r"""A power law function:
    
    :math:`F(E) = A \bigl(\frac{E}{E_{piv}}\bigr)^{\lambda}`,
    
    where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`\lambda` is the power-law index
    * :math:`E_{piv}` is the pivot energy in keV.  Nominally fixed to 100 keV.
    """
    # The total number of parameters
    nparams = 3

    # The list of parameters: (parameter name, units, description)
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('index', '', 'Photon index'),
                  ('Epiv', '', 'Pivot energy')]
    # -----------------------------------------------------------------------
    # The following attributes are used for fitting and parameter estimation

    # The default initialization values for each parameter
    default_values = [0.1, -2.0, 100.0]
    # The maximum absolute fitting step for each parameter 
    delta_abs = [0.1, 0.1, 0.1]
    # The maximum relative fitting step for each parameter 
    delta_rel = [0.01, 0.01, 0.01]
    # The minimum value that each parameter can take
    min_values = [1e-10, -20.0, 0.01]
    # The maximum value that each parameter can take   
    max_values = [np.inf, np.inf, np.inf]
    # If True, then the parameter is a free parameter, otherwise it is fixed  
    free = [True, True, False]

    def eval(self, params, x):
        return params[0] * (x / params[2]) ** params[1]


class Comptonized(Function):
    r"""An exponentially cut-off power law, parametrized by the
    peak of the :math:`\nu F_{\nu}` spectrum:
    
    :math:`F(E) = A e^{-(2+\lambda)E/E_{peak}} \bigl(\frac{E}{E_{piv}}\bigr)^\lambda`
    
    where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`E_{peak}` is the :math:`\nu F_{\nu}` peak in keV
    * :math:`\lambda` is the power-law index
    * :math:`E_{piv}` is the pivot energy in keV.  Nominally fixed to 100 keV.
    """
    nparams = 4
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('Epeak', 'keV', 'SED Peak'),
                  ('index', '', 'Photon index'),
                  ('Epiv', 'keV', 'Pivot energy')]
    default_values = [0.01, 300.0, -0.5, 100.0]
    delta_abs = [0.1, 0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01, 0.01]
    min_values = [1e-10, 0.01, -1.9, 0.01]
    max_values = [np.inf, np.inf, 20.0, np.inf]
    free = [True, True, True, False]

    def eval(self, params, x):
        A, Epeak, index, Epiv = params
        if index == -2.0:
            return np.zeros_like(x)
        logfxn = -x * (2.0 + index) / Epeak + index * np.log(x / Epiv)
        return A * np.exp(logfxn)


class Band(Function):
    r"""The Band GRB function (a type of smoothly broken power law), parametrized 
    by the peak of the :math:`\nu F_{\nu}` spectrum:
    
    :math:`F(E) = A 
    \begin{cases}
    \bigl(\frac{E}{E_{piv}}\bigr)^\alpha e^{-(2+\alpha)E/E_{peak}}, 
    \text{ if } {E < \frac{(\alpha-\beta)E_{peak}}{2+\alpha}}\\
    \bigl(\frac{(\alpha-\beta)E_{peak}}{(2+\alpha)E_{piv}}\bigr)^{\alpha-\beta}
    \exp{\beta-\alpha}\bigl(\frac{E}{E_{piv}}\bigr)^\beta, \text{ otherwise }
    \end{cases}`
    
    where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`E_{peak}` is the :math:`\nu F_{\nu}` peak in keV
    * :math:`\alpha` is the low-energy power-law index
    * :math:`\beta` is the high-energy power-law index
    * :math:`E_{piv}` is the pivot energy in keV.  Nominally fixed to 100 keV.
    """
    nparams = 5
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('Epeak', 'keV', 'SED Peak'),
                  ('alpha', '', 'Low-Energy Photon index'),
                  ('beta', '', 'High-Energy Photon index'),
                  ('Epiv', 'keV', 'Pivot energy')]
    default_values = [0.01, 500.0, -0.5, -2.5, 100.0]
    delta_abs = [0.1, 0.1, 0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01, 0.01, 0.01]
    min_values = [1e-10, 0.01, -1.9, -10.0, 0.01]
    max_values = [np.inf, np.inf, 20.0, -2.0001, np.inf]
    free = [True, True, True, True, False]

    def eval(self, params, x):
        A, Epeak, alpha, beta, Epiv = params
        e0 = Epeak / (2.0 + alpha)
        ebreak = (alpha - beta) * e0
        idx = (x < ebreak)
        logfxn = np.zeros(len(x), dtype=float)
        logfxn[idx] = np.log(A) + alpha * np.log(x[idx] / Epiv) - x[idx] / e0
        dindex = alpha - beta
        logfxn[~idx] = np.log(A) + dindex * np.log(dindex * e0 / Epiv) - dindex + beta * np.log(x[~idx] / Epiv)
        return np.exp(logfxn)


class BandOld(Function):
    r"""The Band GRB function (a type of smoothly broken power law), using the
    original parameterization: 

    :math:`F(E) = A 
    \begin{cases}
    \bigl(\frac{E}{E_{piv}}\bigr)^\alpha e^{-\frac{E}{E_0}}, 
    \text{ if } {E < (\alpha-\beta) E_0}\\
    \bigl(\frac{(\alpha-\beta) E_0}{E_{piv}}\bigr)^{(\alpha-\beta)} 
    e^{(\beta-\alpha)} \bigl(\frac{E}{E_{piv}}\bigr)^\beta, \text{ otherwise }
    \end{cases}`
    
    where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`E_0` is the break energy in keV
    * :math:`\alpha` is the low-energy power-law index
    * :math:`\beta` is the high-energy power-law index
    * :math:`E_{piv}` is the pivot energy in keV.  Nominally fixed to 100 keV.
    
    References:
        `Band, D. et al. 1993, ApJ, 413, 281 
        <http://adsabs.harvard.edu/full/1993ApJ...413..281B>`_
    """
    nparams = 5
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('E0', 'keV', 'Break Energy'),
                  ('alpha', '', 'Low-Energy Photon index'),
                  ('beta', '', 'High-Energy Photon index'),
                  ('Epiv', 'keV', 'Pivot energy')]
    default_values = [0.01, 200.0, -0.5, -2.5, 100.0]
    delta_abs = [0.1, 0.1, 0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01, 0.01, 0.01]
    min_values = [1e-10, 0.01, -1.9, -10.0, 0.01]
    max_values = [np.inf, np.inf, 20.0, -2.0001, np.inf]
    free = [True, True, True, True, False]

    def eval(self, params, x):
        A, e0, alpha, beta, Epiv = params
        ebreak = (alpha - beta) * e0
        idx = (x < ebreak)
        logfxn = np.zeros(len(x), dtype=float)
        logfxn[idx] = np.log(A) + alpha * np.log(x[idx] / Epiv) - x[idx] / e0
        dindex = alpha - beta
        logfxn[~idx] = np.log(A) + dindex * np.log(dindex * e0 / Epiv) - dindex + beta * np.log(x[~idx] / Epiv)
        return np.exp(logfxn)


class BrokenPowerLaw(Function):
    r"""A sharply-broken power law:
    
    :math:`F(E) = A
    \begin{cases}
    \bigl(\frac{E}{E_{piv}}\bigr)^\alpha, & \text{ if } {E \leq E_b}\\
    \bigl(\frac{E_b}{E_{piv}}\bigr)^\alpha \bigl(\frac{E}{E_b}\bigr)^\beta, &
    \text{ otherwise }
    \end{cases}`
    
    where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`E_{b}` is the break energy in keV
    * :math:`\alpha` is the low-energy power-law index
    * :math:`\beta` is the high-energy power-law index
    * :math:`E_{piv}` is the pivot energy in keV.  Nominally fixed to 100 keV.
    """
    nparams = 5
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('Eb', 'keV', 'Break Energy'),
                  ('alpha', '', 'Low-Energy Photon index'),
                  ('beta', '', 'High-Energy Photon index'),
                  ('Epiv', 'keV', 'Pivot energy')]
    default_values = [0.01, 700.0, -1.0, -2.0, 100.0]
    delta_abs = [0.1, 0.1, 0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01, 0.01, 0.01]
    min_values = [1e-10, 0.01, -1.9, -10.0, 0.01]
    max_values = [np.inf, np.inf, 20.0, 20.0, np.inf]
    free = [True, True, True, True, False]

    def eval(self, params, x):
        A, ebreak, alpha, beta, Epiv = params
        mask = (x <= ebreak)
        fxn = np.zeros(len(x), dtype=float)
        fxn[mask] = (x[mask] / Epiv) ** alpha
        fxn[~mask] = (ebreak / Epiv) ** alpha * (x[~mask] / ebreak) ** beta
        return A * fxn


class DoubleBrokenPowerLaw(Function):
    r"""A sharply-broken power law with two breaks:
    
    :math:`F(E) = A
    \begin{cases}
    \bigl(\frac{E}{E_{piv}}\bigr)^\alpha, 
    \text{ if } {E \leq E_{b1}}\\
    \bigl(\frac{E_{b1}}{E_{piv}}\bigr)^\alpha \bigl(\frac{E}{E_b}\bigr)^\beta, 
    \text{ if } {E > E_{b1} \text{ and } E \leq E_{b2}}\\
    \bigl(\frac{E_{b1}}{E_{piv}}\bigr)^\alpha \bigl(\frac{E_{b1}}{E_{b2}}\bigr)^\beta
    \bigl(\frac{E}{E_{b2}}\bigr)^\gamma, 
    \text{ otherwise}
    \end{cases}`        
    
    where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`E_{b1}` is the first break energy in keV
    * :math:`E_{b2}` is the second break energy in keV
    * :math:`\alpha` is the low-energy power-law index
    * :math:`\beta` is the mid-energy power-law index
    * :math:`\gamma` is the high-energy power-law index
    * :math:`E_{piv}` is the pivot energy in keV.  Nominally fixed to 100 keV.
    """
    nparams = 7
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('Eb1', 'keV', 'Lower Break Energy'),
                  ('Eb2', 'keV', 'Upper Break Energy'),
                  ('alpha', '', 'Low-Energy Photon index'),
                  ('beta', '', 'Mid-Energy Photon index'),
                  ('gamma', '', 'High-Energy Photon index'),
                  ('Epiv', 'keV', 'Pivot energy')]
    default_values = [0.01, 200.0, 500, -0.5, -1.0, -2.0, 100.0]
    delta_abs = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    min_values = [1e-10, 0.01, 0.01, -1.9, -10.0, -20.0, 0.01]
    max_values = [np.inf, np.inf, np.inf, 20.0, 20.0, 20.0, np.inf]
    free = [True, True, True, True, True, True, False]

    def eval(self, params, x):
        A, ebreak1, ebreak2, alpha, beta, gamma, Epiv = params
        fxn = np.zeros(len(x), dtype=float)
        mask1 = (x <= ebreak1)
        fxn[mask1] = (x[mask1] / Epiv) ** alpha
        mask2 = (x > ebreak1) & (x <= ebreak2)
        fxn[mask2] = (ebreak1 / Epiv) ** alpha * (x[mask2] / ebreak1) ** beta
        mask3 = (x > ebreak2)
        fxn[mask3] = (ebreak1 / Epiv) ** alpha * (ebreak1 / ebreak2) ** beta * (x[mask3] / ebreak2) ** gamma
        return A * fxn


class SmoothlyBrokenPowerLaw(Function):
    r"""A smoothly-broken power law with a break scale:
    
    :math:`F(E) = A \bigl(\frac{E}{E_{piv}}\bigr)^b 10^{(\beta-\beta_{piv})}, \\
    \text{where }\\
    \beta = m \Delta \ln\bigl(\frac{e^\alpha + e^{-\alpha}}{2}\bigr); \text{    }
    \beta_{piv} = m \Delta \ln\bigl(\frac{e^{\alpha_{piv}} + e^{-\alpha_{piv}}}{2}\bigr); \\
    \alpha = \frac{\log(E/E_b)}{\Delta}; \text{    }
    \alpha_{piv} = \frac{\log(E_{piv}/E_b)}{\Delta}; \\
    m = \frac{\lambda_2 - \lambda_1}{2}; \text{    }
    b = \frac{\lambda_1 + \lambda_2}{2};`
    
    and where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`E_{b}` is the break energy in keV
    * :math:`\Delta` is the break scale in decades of energy
    * :math:`\lambda_1` is the low-energy power-law index
    * :math:`\lambda_2` is the high-energy power-law index
    * :math:`E_{piv}` is the pivot energy in keV.  Nominally fixed to 100 keV.
    """
    nparams = 6
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('Eb', 'keV', 'Break Energy'),
                  ('Delta', 'dex', 'Break Scale'),
                  ('lambda1', '', 'Low-Energy Photon index'),
                  ('lambda2', '', 'High-Energy Photon index'),
                  ('Epiv', 'keV', 'Pivot energy')]
    default_values = [0.01, 700.0, 0.3, -1.0, -2.0, 100.0]
    delta_abs = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    min_values = [1e-10, 0.01, 0.01, -1.9, -10.0, 0.01]
    max_values = [np.inf, np.inf, 10.0, 20.0, 20.0, np.inf]
    free = [True, True, True, True, True, False]

    def eval(self, params, x):
        A, ebreak, delta, lam1, lam2, Epiv = params
        m = (lam2 - lam1) / 2.0
        b = (lam1 + lam2) / 2.0
        apiv = np.log10(Epiv / ebreak) / delta
        a = np.log10(x / ebreak) / delta
        bpiv = m * delta * np.log((np.exp(apiv) + np.exp(-apiv)) / 2.0)
        beta = m * delta * np.log((np.exp(a) + np.exp(-a)) / 2.0)
        fxn = A * (x / Epiv) ** b * 10.0 ** (beta - bpiv)
        return fxn


class LogNormal(Function):
    r"""A Log-Normal function:
        
    :math:`F(E) = \frac{A}{\sqrt{2\pi} \sigma E} 
    e^{-\bigl(\frac{\ln E - \mu}{2\sigma}\bigr)^2}`
    
    where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`\mu` is the natural log of energy
    * :math:`\sigma` is the logarithmic standard deviation
    """
    nparams = 3
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('mu', 'Ln(keV)', 'Ln(Energy)'),
                  ('sigma', 'Ln(keV)', 'Log Standard Deviation')]
    default_values = [0.01, 5.0, 1.0]
    delta_abs = [0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01]
    min_values = [1e-10, -1.0, 0.0]
    max_values = [np.inf, 10.0, np.inf]
    free = [True, True, True]

    def eval(self, params, x):
        A, mu, sig = params
        fxn = A / (np.sqrt(2.0 * np.pi) * sig * x) * np.exp(-0.5 * ((np.log(x) - mu) / sig) ** 2)
        return fxn


class GaussianLog(Function):
    r"""A Gaussian with :math:`log_{10}(E)`:
        
    :math:`F(E) = \frac{A}{\sqrt{2\pi} \sigma} 
    e^{-\bigl(\frac{\log E-\log E_{cen} }{2\sigma}\bigr)^2}, \\
    \text{ where } \\
    \sigma = \frac{W}{2.35482}`
    
    and where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`E_{cen}` is the centroid energy in keV
    * :math:`W` is the :math:`\log_{10}` FWHM at :math:`E_{cen}` in decades
      of energy
    """
    nparams = 3
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('Ecen', 'keV', 'Centroid'),
                  ('W', 'dex', 'Log(FWHM)')]
    default_values = [0.01, 100.0, 1.0]
    delta_abs = [0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01]
    min_values = [1e-10, 0.1, 0.1]
    max_values = [np.inf, np.inf, np.inf]
    free = [True, True, True]

    def eval(self, params, x):
        A, Ecen, w = params
        sig = w / 2.35482
        fxn = A / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-0.5 * ((np.log10(x) - np.log10(Ecen)) / sig) ** 2)
        return fxn


class GaussianLogVaryingFWHM(Function):
    r"""A Gaussian with :math:`log_{10}(E)` and linearly-varying FWHM:
        
    :math:`F(E) = \frac{A}{\sqrt{2\pi} \sigma} 
    e^{-\bigl(\frac{\log E-\log E_{cen} }{2\sigma}\bigr)^2}, \\
    \text{ where } \\
    \sigma = \frac{W + s (\log E - \log E_{cen})}{2.35482}`
    
    and where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`E_{cen}` is the centroid energy in keV
    * :math:`W` is the :math:`\log_{10}` FWHM at :math:`E_{cen}` in decades
      of energy
    * :math:`s` is the slope of :math:`W` with respect to :math:`\log_{10}(E)`
      in decades per  :math:`\log_{10}(E)`
    """
    nparams = 4
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('Ecen', 'keV', 'Centroid'),
                  ('W', 'dex', 'Log(FWHM)'),
                  ('s', 'dex/log(keV)', 'Slope of W')]
    default_values = [0.01, 100.0, 1.0, 0.0]
    delta_abs = [0.1, 0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01, 0.1]
    min_values = [1e-10, 0.1, 0.1, -1.0]
    max_values = [np.inf, np.inf, np.inf, np.inf]
    free = [True, True, True, True]

    def eval(self, params, x):
        A, Ecen, w, s = params
        sig = (w + s * (np.log10(x) - np.log10(Ecen))) / 2.35482
        fxn = A / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-0.5 * ((np.log10(x) - np.log10(Ecen)) / sig) ** 2)
        return fxn


class SunyaevTitarchuk(Function):
    r"""Sunyaev-Titarchuk Comptonization spectra:
        
    :math:`F(E) = A E^2 e^{-E} \int_0^\infty t^{n-1} 
    e^{-t \bigl(1+\frac{t}{E}\bigr)^{n+3}} \text{d}t, \\
    \text{ where } \\
    n = -\frac{3}{2} + \sqrt{\gamma + \frac{9}{4}}; \\
    \gamma = \frac{511 \pi^2}{G kT \bigl(\tau + \frac{2}{3} \bigr)^2}`
    
    and where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`kT` is the electron energy in keV
    * :math:`\tau` is the optical depth
    * :math:`G` is the geometry factor, fixed at 3 for a spherical cloud and
      at 12 for a disk of electrons
      
    References:
        `Sunyaev, R. A. and Titarchuk, L. G. 1980, A&A, 86, 121 
        <http://adsabs.harvard.edu/full/1980A%26A....86..121S>`_
    """
    nparams = 4
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('kT', 'keV', 'Electron energy'),
                  ('tau', '', 'Optical depth'),
                  ('G', '', 'Geometry factor')]
    default_values = [0.01, 30.0, 10.0, 3.0]
    delta_abs = [0.1, 0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01, 0.01]
    min_values = [1e-10, 0.1, 0.0, 3]
    max_values = [np.inf, np.inf, np.inf, 12]
    free = [True, True, True, False]

    def eval(self, params, x):
        A, kT, tau, G = params
        gamma = 511.0 * np.pi ** 2 / (G * kT * (tau + 2. / 3.) ** 2)
        n = -1.5 + np.sqrt(gamma + 2.25)

        def integrand(t_v, n_v, x_v):
            return t_v ** (n_v - 1.0) * np.exp(-t_v * (1.0 + t_v / x_v) ** (n_v + 3.))

        integral = [quad(integrand, 0, np.inf, args=(n, x1))[0] for x1 in x]
        fxn = A * x ** 2 * np.exp(-x) * np.array(integral)
        return fxn


class OTTB(Function):
    r"""Optically-Thin Thermal Bremsstrahlung:
        
    :math:`F(E) = A e^{-\frac{E}{kT}} e^{\frac{E_{piv}}{kT}} 
    \biggl(\frac{E_{piv}}{E} \biggr)`    
    
    where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`kT` is the electron energy in keV
    * :math:`E_{piv}` is the pivot energy in keV
    """
    nparams = 3
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('kT', 'keV', 'Electron energy'),
                  ('Epiv', 'keV', 'Pivot energy')]
    default_values = [0.01, 30.0, 100.0]
    delta_abs = [0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01]
    min_values = [1e-10, 0.1, 0.0]
    max_values = [np.inf, np.inf, np.inf]
    free = [True, True, False]

    def eval(self, params, x):
        A, kT, Epiv = params
        fxn = A * np.exp(-x / kT) * np.exp(Epiv / kT) * (Epiv / x)
        return fxn


class BlackBody(Function):
    r"""Black Body spectrum:
        
    :math:`F(E) = A \frac{E^2}{e^{E/kT}-1}`    
    
    where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`kT` is the temperature in keV
    """
    nparams = 2
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('kT', 'keV', 'Temperature')]
    default_values = [0.01, 30.0]
    delta_abs = [0.1, 0.1]
    delta_rel = [0.01, 0.01]
    min_values = [1e-10, 0.1]
    max_values = [np.inf, np.inf]
    free = [True, True]

    def eval(self, params, x):
        A, kT = params
        fxn = A * x ** 2 / (np.exp(x / kT) - 1.0)
        return fxn


class YangSoongPulsar(Function):
    r"""Yang Soong's pulsar spectral form
        
    :math:`F(E) = A C (1 - E_W G), \\
    \text{ where} \\
    C = \begin{cases}
    \bigl(\frac{E}{E_{piv}}\bigr)^{-\alpha}, & \text{ if } {E \leq E_b} \\
    \frac{M}{E} e^{-E/E_{piv}}, & \text{ otherwise }
    \end{cases}; \\
    M = e^{E_b/E_F} E_b \bigl(\frac{E_b}{E_{piv}}\bigr)^{-\alpha}; \\
    G = \frac{0.94}{W} e^{-2.76 (\frac{E-E_{cen}}{W})^2};`
    
    and where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`\alpha` is the power-law index
    * :math:`E_b` is the break energy in keV
    * :math:`E_F` is the e-folding energy in keV
    * :math:`E_W` is the equivalent line width in keV
    * :math:`E_{cen}` is the line centroid in keV
    * :math:`W` is the FWHM of the line, in keV
    * :math:`E_{piv}` is the pivot energy in keV
    
    References:
        `Soong, Y., et al. 1990, ApJ, 348, 641 
        <http://adsabs.harvard.edu/full/1990ApJ...348..641S>`_
    """
    nparams = 8
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('alpha', '', 'Index'),
                  ('Eb', 'keV', 'Break Energy'),
                  ('EF', 'keV', 'E-folding Energy'),
                  ('EW', 'keV', 'Line Width'),
                  ('Ecen', 'keV', 'Line Centroid'),
                  ('W', 'keV', 'Line FWHM'),
                  ('Epiv', 'keV', 'Pivot Energy')]
    default_values = [0.01, -1.0, 200.0, 300.0, 50.0, 200.0, 50.0, 100.0]
    delta_abs = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    min_values = [1e-10, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    max_values = [np.inf, 10.0, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]
    free = [True, True, True, True, True, True, True, False]

    def eval(self, params, x):
        A, alpha, Eb, EF, EW, Ecen, W, Epiv = params
        M = np.exp(Eb / EF) * Eb * (Eb / Epiv) ** (-alpha)
        G = 0.94 / W * np.exp(-2.76 * ((x - Ecen) / W) ** 2)

        C = np.zeros(len(x), dtype=float)
        mask = (x <= Eb)
        C[mask] = A * (x[mask] / Epiv) ** (-alpha)
        C[~mask] = A * M * np.exp(-x[~mask] / Epiv) / x[~mask]
        fxn = C * (1.0 - EW * G)
        return fxn


class TanakaPulsar(Function):
    r"""Tanaka's pulsar model
        
    :math:`F(E) = A \frac{C}{C_{piv}} e^{L_{piv}-L}, \\
    \text{ where} \\
    C = A E^{-\alpha} e^{-E/kT}; \text{    }
    C_{piv} E_{piv}^{-\alpha} e^{-E_{piv}/kT}; \\
    L = \sum_{i=1}^{N_L} \frac{\tau_i (i W)^2 [E/(i E_L)]^2}{(E-i E_L)^2 + (i W)^2}; \\
    L_{piv} = \sum_{i=1}^{N_L} \frac{\tau_i (i W)^2 [E_{piv}/(i E_L)]^2} 
    {(E_{piv}-i E_L)^2 + (i W)^2}; \\`
    
    and where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`\alpha` is the power-law index
    * :math:`kT` is the temperature in keV
    * :math:`\tau_1` is the optical depth of the first line
    * :math:`\tau_2` is the optical depth of the second line
    * :math:`N_L` is the number of lines; fixed at either 0, 1, or 2
    * :math:`E_L` is the line centroid in keV of the first line
    * :math:`W` is the line FWHM in keV
    * :math:`E_{piv}` is the pivot energy in keV
    
    
    References:
        `Grove, J. E., et al. 1995, ApJ, 438, L25 
        <https://apps.dtic.mil/docs/citations/ADA461878>`_
    """
    nparams = 9
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('alpha', '', 'Index'),
                  ('kT', 'keV', 'Temperature'),
                  ('tau1', '', 'Optical Depth of 1st line'),
                  ('tau2', '', 'Optical Depth of 2nd line'),
                  ('NL', '', 'Number of lines'),
                  ('EL', 'keV', '1st Line Centroid'),
                  ('W', 'keV', '1st Line FWHM'),
                  ('Epiv', 'keV', 'Pivot Energy')]
    default_values = [0.01, -1.0, 200.0, 1.0, 1.0, 1, 100.0, 50.0, 100.0]
    delta_abs = [0.1, 0.1, 0.1, 0.1, 0.1, 0.0, 0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01, 0.01, 0.01, 0.00, 0.01, 0.01, 0.01]
    min_values = [1e-10, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    max_values = [np.inf, 10.0, np.inf, np.inf, np.inf, 2, np.inf, np.inf, np.inf]
    free = [True, True, True, True, True, False, True, True, False]

    def eval(self, params, x):
        A, alpha, kT, tau1, tau2, NL, EL, W, Epiv = params
        C = A * x ** (-alpha) * np.exp(-x / kT)
        Cpiv = Epiv ** (-alpha) * np.exp(-Epiv / kT)

        if NL > 0:
            taus = np.array((tau1, tau2))
            idx = np.arange(2) + 1

            L = ((taus * (idx[np.newaxis, :] * W) ** 2 * (x[:, np.newaxis] / (idx[np.newaxis, :] * EL)) ** 2)
                 / ((x[:, np.newaxis] - idx[np.newaxis, :] * EL) ** 2 + (idx[np.newaxis, :] * W) ** 2))
            L = L[:, :NL].sum(axis=1)
            Lpiv = (taus * (idx * W) ** 2 * (Epiv / (idx * EL)) ** 2) / ((Epiv - idx * EL) ** 2 + (idx * W) ** 2)
            Lpiv = Lpiv[:NL].sum()
        else:
            L = Lpiv = 0.0

        fxn = A * C / Cpiv * np.exp(Lpiv - L)
        return fxn


class OTTS(Function):
    r"""Optically-Thin Thermal Synchrotron:
        
    :math:`F(E) = A e^{(E/E_C)^{1/3}}`
    
    where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`E_C` is the energy scale in keV
    
    References:
        `Liang, E. P., et al., 1983, ApJ, 271, 776 
        <http://adsabs.harvard.edu/full/1983ApJ...271..766L>`_
    """
    nparams = 2
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('EC', 'keV', 'Energy Scale')]
    default_values = [0.01, 100.0]
    delta_abs = [0.1, 0.1]
    delta_rel = [0.01, 0.01]
    min_values = [1e-10, 0]
    max_values = [np.inf, np.inf]
    free = [True, True]

    def eval(self, params, x):
        A, Ec = params
        fxn = A * np.exp((x / Ec) ** (1.0 / 3.0))
        return fxn


class GaussLine(Function):
    r"""A Gaussian Line:
        
    :math:`F(E) = A \frac{\text{erfc}(a_l)-\text{erfc}(a_r)}{2(E_r-E_l)}, \\
    \text{ where } \\
    a_l = \frac{E_l - E_{cen}}{\sqrt{2}\sigma}; \text{    }
    a_r = \frac{E_r - E_{cen}}{\sqrt{2}\sigma}; \\
    \sigma = W/2.35482`
    
    and where 
    
    * :math:`A` is the amplitude in ph/s/cm^2/keV
    * :math:`E_{cen}` is the centroid energy in keV
    * :math:`W` is the FHWM in keV
    """
    nparams = 3
    param_list = [('A', 'ph/s/cm^2/keV', 'Amplitude'),
                  ('Ecen', 'keV', 'Centroid Energy'),
                  ('W', 'keV', 'FWHM')]
    default_values = [0.01, 100.0, 8.0]
    delta_abs = [0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01]
    min_values = [1e-10, 0.0, 0.0]
    max_values = [np.inf, np.inf, np.inf]
    free = [True, True, True]

    def eval(self, params, x):
        A, Ecen, W = params

        # calculate edges; need to add an ending edge
        edges = np.asarray(x)
        dE = np.log10(x[-1]) - np.log10(x[-2])
        edges = np.append(edges, [x[-1] * 10.0 ** dE])
        el = edges[:-1]
        er = edges[1:]

        sig = W / 2.35482
        al = (el - Ecen) / (np.sqrt(2.0) * sig)
        ar = (er - Ecen) / (np.sqrt(2.0) * sig)
        fxn = A * (erfc(al) - erfc(ar)) / (2.0 * (er - el))
        return fxn


class LowEnergyCutoff(Function):
    r"""A Low-Energy cutoff (Multiplicative Component):
        
    :math:`F(E) = \begin{cases}
    (E/E_{cut})^{E_{cut}/E_F} e^{(E_{cut}-E)/E_F} & \text{ if } {E \leq E_{cut}} \\
    1 & \text{ otherwise }
    \end{cases}`
    
    where 
    
    * :math:`E_{cut}` is the cutoff energy in keV
    * :math:`E_{cen}` is the folding energy in keV
    """
    nparams = 2
    param_list = [('Ecut', 'keV', 'Cutoff Energy'),
                  ('EF', 'keV', 'Folding Energy')]
    default_values = [100.0, 10.0]
    delta_abs = [0.1, 0.1]
    delta_rel = [0.01, 0.01]
    min_values = [0.0, 0.0]
    max_values = [np.inf, np.inf]
    free = [True, True]

    def eval(self, params, x):
        Ecut, EF = params
        mask = (x <= Ecut)
        fxn = np.ones(len(x), dtype=float)
        fxn[mask] = (x[mask] / Ecut) ** (Ecut / EF) * np.exp((Ecut - x[mask]) / EF)
        return fxn


class HighEnergyCutoff(Function):
    r"""A High-Energy cutoff (Multiplicative Component):
        
    :math:`F(E) = \begin{cases}    
    (E/E_{cut})^{E_{cut}/E_F} e^{(E_{cut}-E)/E_F} & \text{ if } {E > E_{cut}} \\
    1 & \text{ otherwise }
    \end{cases}`
    
    where 
    
    * :math:`E_{cut}` is the cutoff energy in keV
    * :math:`E_{cen}` is the folding energy in keV
    """
    nparams = 2
    param_list = [('Ecut', 'keV', 'Cutoff Energy'),
                  ('EF', 'keV', 'Folding Energy')]
    default_values = [1000.0, 100.0]
    delta_abs = [0.1, 0.1]
    delta_rel = [0.01, 0.01]
    min_values = [0.0, 0.0]
    max_values = [np.inf, np.inf]
    free = [True, True]

    def eval(self, params, x):
        Ecut, EF = params
        mask = (x > Ecut)
        fxn = np.ones(len(x), dtype=float)
        fxn[mask] = (x[mask] / Ecut) ** (Ecut / EF) * np.exp((Ecut - x[mask]) / EF)
        return fxn


class PowerLawMult(Function):
    r"""A Power Law (Multiplicative Component):
        
    :math:`F(E) = \bigl(\frac{E}{E_{piv}} \bigr)^\lambda`
    
    where 
    
    * :math:`\lambda` is the index
    * :math:`E_{piv}` is the pivot energy in keV
    """
    nparams = 2
    param_list = [('lambda', '', 'Index'),
                  ('Epiv', 'keV', 'Pivot Energy')]
    default_values = [-2.0, 100.0]
    delta_abs = [0.1, 0.1]
    delta_rel = [0.01, 0.01]
    min_values = [-10.0, 0.0]
    max_values = [10.0, np.inf]
    free = [True, False]

    def eval(self, params, x):
        lam, Epiv = params
        fxn = (x / Epiv) ** lam
        return fxn


class GaussLineMult(Function):
    r"""A Gaussian Line (Multiplicative Component):
        
    :math:`F(E) = e^{I G}, \\
    \text{ where } \\
    G = \frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{(E-E_{cen})^2}{2\sigma^2}}; \\
    \sigma = W/2.35482`
    
    and where 
    
    * :math:`I` is the intensity
    * :math:`E_{cen}` is the centroid energy in keV
    * :math:`W` is the FWHM in keV    
    """
    nparams = 3
    param_list = [('I', '', 'Intensity'),
                  ('Ecen', 'keV', 'Centroid Energy'),
                  ('W', 'keV', 'FWHM')]
    default_values = [1.0, 10.0, 8.0]
    delta_abs = [0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01]
    min_values = [0.0, 0.0, 0.0]
    max_values = [np.inf, np.inf, np.inf]
    free = [True, True, True]

    def eval(self, params, x):
        I, Ecen, W = params
        sig = W / 2.35482
        G = 1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-(x - Ecen) ** 2 / (2.0 * sig ** 2))
        fxn = np.exp(I * G)
        return fxn


class LorentzLineMult(Function):
    r"""A Lorentzian Line (Multiplicative Component):
        
    :math:`F(E) = e^{\frac{I}{W^2 + (E-E_{cen})^2}}`
    
    where 
    
    * :math:`I` is the intensity
    * :math:`E_{cen}` is the centroid energy in keV
    * :math:`W` is the FWHM in keV    
    """
    nparams = 3
    param_list = [('I', '', 'Intensity'),
                  ('Ecen', 'keV', 'Centroid Energy'),
                  ('W', 'keV', 'FWHM')]
    default_values = [1.0, 10.0, 8.0]
    delta_abs = [0.1, 0.1, 0.1]
    delta_rel = [0.01, 0.01, 0.01]
    min_values = [0.0, 0.0, 0.0]
    max_values = [np.inf, np.inf, np.inf]
    free = [True, True, True]

    def eval(self, params, x):
        I, Ecen, W = params
        fxn = np.exp(I / (W ** 2 + (x - Ecen) ** 2))
        return fxn
