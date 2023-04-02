#
#     Authors: William Cleveland (USRA),
#              Adam Goldstein (USRA) and
#              Daniel Kocevski (NASA)
#
#     Portions of the code are Copyright 2020 William Cleveland and
#     Adam Goldstein, Universities Space Research Association
#     All rights reserved.
#
#     Written for the Fermi Gamma-ray Burst Monitor (Fermi-GBM)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import unittest
import numpy as np
from gdt.core.spectra.parameters import *
from gdt.core.data_primitives import Parameter


class TestPhotonFlux(unittest.TestCase):

    def setUp(self):
        self.pflux = PhotonFlux(2.0, (0.1, 0.2), (50.0, 300.0))

    def test_name(self):
        self.assertEqual(self.pflux.name, 'Photon Flux')

    def test_support(self):
        self.assertTupleEqual(self.pflux.support, (0.0, np.inf))

    def test_uncertainty(self):
        self.assertTupleEqual(self.pflux.uncertainty, (0.1, 0.2))

    def test_units(self):
        self.assertEqual(self.pflux.units, 'ph/cm^2/s')

    def test_value(self):
        self.assertEqual(self.pflux.value, 2.0)

    def test_energy_range(self):
        self.assertTupleEqual(self.pflux.energy_range, (50.0, 300.0))


class TestPhotonFluence(unittest.TestCase):

    def setUp(self):
        self.pfluence = PhotonFluence(20.0, (0.4, 0.8), (50.0, 300.0))

    def test_name(self):
        self.assertEqual(self.pfluence.name, 'Photon Fluence')

    def test_support(self):
        self.assertTupleEqual(self.pfluence.support, (0.0, np.inf))

    def test_uncertainty(self):
        self.assertTupleEqual(self.pfluence.uncertainty, (0.4, 0.8))

    def test_units(self):
        self.assertEqual(self.pfluence.units, 'ph/cm^2')

    def test_value(self):
        self.assertEqual(self.pfluence.value, 20.0)

    def test_energy_range(self):
        self.assertTupleEqual(self.pfluence.energy_range, (50.0, 300.0))


class TestEnergyFlux(unittest.TestCase):

    def setUp(self):
        self.eflux = EnergyFlux(2e-7, (1e-8, 2e-8), (50.0, 300.0))

    def test_name(self):
        self.assertEqual(self.eflux.name, 'Energy Flux')

    def test_support(self):
        self.assertTupleEqual(self.eflux.support, (0.0, np.inf))

    def test_uncertainty(self):
        self.assertTupleEqual(self.eflux.uncertainty, (1e-8, 2e-8))

    def test_units(self):
        self.assertEqual(self.eflux.units, 'erg/cm^2/s')

    def test_value(self):
        self.assertEqual(self.eflux.value, 2e-7)

    def test_energy_range(self):
        self.assertTupleEqual(self.eflux.energy_range, (50.0, 300.0))


class TestEnergyFluence(unittest.TestCase):

    def setUp(self):
        self.efluence = EnergyFluence(2e-6, (4e-8, 8e-8), (50.0, 300.0))

    def test_name(self):
        self.assertEqual(self.efluence.name, 'Energy Fluence')

    def test_support(self):
        self.assertTupleEqual(self.efluence.support, (0.0, np.inf))

    def test_uncertainty(self):
        self.assertTupleEqual(self.efluence.uncertainty, (4e-8, 8e-8))

    def test_units(self):
        self.assertEqual(self.efluence.units, 'erg/cm^2')

    def test_value(self):
        self.assertEqual(self.efluence.value, 2e-6)

    def test_energy_range(self):
        self.assertTupleEqual(self.efluence.energy_range, (50.0, 300.0))


class TestModelFit(unittest.TestCase):

    def setUp(self):
        amp = Parameter(0.01, 2e-3, name='A', units='ph/cm^2/s/keV')
        index = Parameter(-1.7, 0.3, name='index')
        self.params = [amp, index]
        self.model = ModelFit.from_data('PowerLaw', (0.0, 10.0),
                                        parameters=self.params)

    def test_covariance(self):
        self.model.covariance = np.array([[0.2, 1.0], [1.0, 0.2]])
        self.assertListEqual(self.model.covariance.flatten().tolist(),
                             [0.2, 1.0, 1.0, 0.2])

        # wrong type           
        with self.assertRaises(TypeError):
            self.model.covariance = True

        # wrong number of dimensions
        with self.assertRaises(ValueError):
            self.model.covariance = np.array([0.2, 1.0])

        # wrong shape
        with self.assertRaises(ValueError):
            self.model.covariance = np.array([[0.2, 1.0, 0.0], [1.0, 0.2, 0.0]])

    def test_dof(self):
        self.model.dof = 22
        self.assertEqual(self.model.dof, 22)

        # wrong type
        with self.assertRaises(TypeError):
            self.model.dof = 'oops'

    def test_energy_fluence(self):
        self.model.energy_fluence = EnergyFluence(2e-6, 1e-7, (50.0, 300.0))
        self.assertEqual(self.model.energy_fluence.value, 2e-6)
        self.assertTupleEqual(self.model.energy_fluence.uncertainty, (1e-7, 1e-7))

        # wrong type
        with self.assertRaises(TypeError):
            self.model.energy_fluence = (2e-6, 1e-7)

    def test_energy_flux(self):
        self.model.energy_flux = EnergyFlux(2e-7, 1e-8, (50.0, 300.0))
        self.assertEqual(self.model.energy_flux.value, 2e-7)
        self.assertTupleEqual(self.model.energy_flux.uncertainty, (1e-8, 1e-8))

        # wrong type
        with self.assertRaises(TypeError):
            self.model.energy_flux = (2e-7, 1e-8)

    def test_flux_energy_range(self):
        self.model.flux_energy_range = (50.0, 300.0)
        self.assertTupleEqual(self.model.flux_energy_range, (50.0, 300.0))

        # wrong type
        with self.assertRaises(TypeError):
            self.model.flux_energy_range = 50.0

        # wrong length
        with self.assertRaises(ValueError):
            self.model.flux_energy_range = (50.0, 300.0, 500.0)

    def test_name(self):
        self.assertEqual(self.model.name, 'PowerLaw')

    def test_parameters(self):
        self.model.parameters = self.params
        self.assertEqual(self.model.parameters[0].value, 0.01)
        self.assertEqual(self.model.parameters[1].value, -1.7)

        # must be a list
        with self.assertRaises(TypeError):
            self.model.parameters = self.params[0]

        # wrong type
        with self.assertRaises(TypeError):
            self.model.parameters = (0.01, -1.7)

    def test_photon_fluence(self):
        self.model.photon_fluence = PhotonFluence(20.0, 0.1, (50.0, 300.0))
        self.assertEqual(self.model.photon_fluence.value, 20.0)
        self.assertTupleEqual(self.model.photon_fluence.uncertainty, (0.1, 0.1))

        # wrong type
        with self.assertRaises(TypeError):
            self.model.photon_fluence = (2.0, 0.1)

    def test_photon_flux(self):
        self.model.photon_flux = PhotonFlux(2.0, 0.1, (50.0, 300.0))
        self.assertEqual(self.model.photon_flux.value, 2.0)
        self.assertTupleEqual(self.model.photon_flux.uncertainty, (0.1, 0.1))

        # wrong type
        with self.assertRaises(TypeError):
            self.model.photon_flux = (2.0, 0.1)

    def test_stat_name(self):
        self.model.stat_name = 'C-stat'
        self.assertEqual(self.model.stat_name, 'C-stat')

    def test_stat_value(self):
        self.model.stat_value = 23.45
        self.assertEqual(self.model.stat_value, 23.45)

        # wrong type
        with self.assertRaises(TypeError):
            self.model.stat_value = 'oops'

    def test_time_range(self):
        self.assertTupleEqual(self.model.time_range, (0.0, 10.0))

    def test_parameter_list(self):
        self.assertListEqual(self.model.parameter_list(), ['A', 'index'])

    def test_init_errors(self):
        # time range wrong type
        with self.assertRaises(TypeError):
            ModelFit.from_data('PowerLaw', 0.0)

        # time range wrong length
        with self.assertRaises(ValueError):
            ModelFit.from_data('PowerLaw', (0.0, 10.0, 20.0))

        # try to initialize an invalid attribute
        with self.assertRaises(AttributeError):
            ModelFit.from_data('PowerLaw', (0.0, 10.0), foo='bar')

    def test_repr(self):
        expected = f'<ModelFit: {self.model.name}>'
        self.assertEqual(expected, repr(self.model))

    def test_str(self):
        param_str = '\n   '.join([str(param) for param in self.model.parameters])
        expected = f'{self.model.name}\n   {param_str}'
        self.assertEqual(expected, str(self.model))


class TestDetectorData(unittest.TestCase):

    def setUp(self):
        self.data = DetectorData.from_data('Fermi, GBM', 'NAI_00', 'TTE', 128,
                                           active=True)

    def test_active(self):
        self.data.active = False
        self.assertFalse(self.data.active)

    def test_channel_mask(self):
        self.data.channel_mask = np.array([1, 0, 1, 0], dtype=bool)
        self.assertListEqual(self.data.channel_mask.tolist(), [1, 0, 1, 0])

        # wrong type
        with self.assertRaises(TypeError):
            self.data.channel_mask = [1, 0, 1, 0]

    def test_channel_range(self):
        self.data.channel_range = (4, 126)
        self.assertTupleEqual(self.data.channel_range, (4, 126))

        # wrong type
        with self.assertRaises(TypeError):
            self.data.channel_range = 4

        # wrong length
        with self.assertRaises(ValueError):
            self.data.channel_range = (4, 5, 126)

        # wrong order        
        with self.assertRaises(ValueError):
            self.data.channel_range = (126, 4)

    def test_datatype(self):
        self.assertEqual(self.data.datatype, 'TTE')

    def test_detector(self):
        self.assertEqual(self.data.detector, 'NAI_00')

    def test_energy_edges(self):
        self.data.energy_edges = np.array([5.0, 10.0, 15.0, 20.0])
        self.assertListEqual(self.data.energy_edges.tolist(),
                             [5.0, 10.0, 15.0, 20.0])

        # wrong type
        with self.assertRaises(TypeError):
            self.data.energy_edges = [5.0, 10.0, 15.0, 20.0]

    def test_energy_range(self):
        self.data.energy_range = (8.0, 900.0)
        self.assertTupleEqual(self.data.energy_range, (8.0, 900.0))

        # wrong type
        with self.assertRaises(TypeError):
            self.data.energy_range = 8.0

        # wrong length
        with self.assertRaises(ValueError):
            self.data.energy_range = (8.0, 50.0, 900.0)

        # wrong order        
        with self.assertRaises(ValueError):
            self.data.energy_range = (900.0, 8.0)

    def test_filename(self):
        self.data.filename = 'myfile.dat'
        self.assertEqual(self.data.filename, 'myfile.dat')

    def test_instrument(self):
        self.assertEqual(self.data.instrument, 'Fermi, GBM')

    def test_num_chans(self):
        self.assertEqual(self.data.num_chans, 128)

    def test_photon_counts(self):
        self.data.photon_counts = np.array([10.5, 11.0, 3.0, 0.0])
        self.assertListEqual(self.data.photon_counts.tolist(), [10.5, 11.0, 3.0, 0.0])

        # wrong type
        with self.assertRaises(TypeError):
            self.data.photon_counts = [10.5, 11.0, 3.0, 0.0]

    def test_photon_errors(self):
        self.data.photon_errors = np.array([1.1, 1.1, 0.3, 0.0])
        self.assertListEqual(self.data.photon_errors.tolist(), [1.1, 1.1, 0.3, 0.0])

        # wrong type
        with self.assertRaises(TypeError):
            self.data.photon_errors = [1.1, 1.1, 0.3, 0.0]

    def test_photon_model(self):
        self.data.photon_model = np.array([10.5, 11.0, 3.0, 0.0])
        self.assertListEqual(self.data.photon_model.tolist(), [10.5, 11.0, 3.0, 0.0])

        # wrong type
        with self.assertRaises(TypeError):
            self.data.photon_model = [10.5, 11.0, 3.0, 0.0]

    def test_response(self):
        self.data.response = 'response.rsp'
        self.assertEqual(self.data.response, 'response.rsp')

    def test_time_range(self):
        self.data.time_range = (0.0, 10.0)
        self.assertTupleEqual(self.data.time_range, (0.0, 10.0))

        # wrong type
        with self.assertRaises(TypeError):
            self.data.time_range = 0.0

        # wrong length
        with self.assertRaises(ValueError):
            self.data.time_range = (0.0, 5.0, 10.0)

        # wrong order        
        with self.assertRaises(ValueError):
            self.data.time_range = (10.0, 0.0)

    def test_repr(self):
        expected = (f'<DetectorData: {self.data.detector}; {self.data.datatype};\n'
                    f' time range: {self.data.time_range} s;\n'
                    f' energy range: {self.data.energy_range} keV>')
        self.assertEqual(expected, repr(self.data))

    def test_errors(self):
        # try to initialize an invalid attribute
        with self.assertRaises(AttributeError):
            DetectorData.from_data('Fermi, GBM', 'NAI_00', 'TTE', 128, foo='bar')
