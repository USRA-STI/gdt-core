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
import os
import unittest
import numpy as np
import numpy.testing as npt

from gdt.core.spectra.fitting import *
from gdt.core.spectra.functions import PowerLaw, BlackBody

from gdt.core.background.primitives import BackgroundSpectrum
from gdt.core.data_primitives import EnergyBins, Gti, ResponseMatrix
from gdt.core.pha import Pha, Bak
from gdt.core.response import Rsp


# power law (0.05, -1.3)
def make_first_pha():
    counts = [235, 276, 212, 24]
    emin = [4.6, 27.3, 102., 538.]
    emax = [27.3, 102., 538., 2000.]
    exposure = 0.256
    data = EnergyBins(counts, emin, emax, exposure)
    gti = Gti.from_list([(0.0, 0.256)])
    pha = Pha.from_data(data, gti=gti)
    return pha


def make_second_pha():
    counts = [255, 283, 228, 22]
    emin = [4.6, 27.3, 102., 538.]
    emax = [27.3, 102., 538., 2000.]
    exposure = 0.256
    data = EnergyBins(counts, emin, emax, exposure)
    gti = Gti.from_list([(0.0, 0.256)])
    pha = Pha.from_data(data, gti=gti)
    return pha


def make_first_bak():
    rates = [37.4443041, 53.72757004, 16.43976248, 31.63717691]
    uncert = [1.896, 2.889, 0.919, 1.66]
    emin = [4.6, 27.3, 102., 538.]
    emax = [27.3, 102., 538., 2000.]
    exposure = 0.256
    data = BackgroundSpectrum(rates, uncert, emin, emax, exposure)
    gti = Gti.from_list([(0.0, 0.256)])
    bak = Bak.from_data(data, gti=gti)
    return bak


def make_second_bak():
    rates = [40.61084827, 55.21141574, 19.38057642, 33.97849876]
    uncert = [1.896, 2.889, 0.919, 1.66]
    emin = [4.6, 27.3, 102., 538.]
    emax = [27.3, 102., 538., 2000.]
    exposure = 0.256
    data = BackgroundSpectrum(rates, uncert, emin, emax, exposure)
    gti = Gti.from_list([(0.0, 0.256)])
    bak = Bak.from_data(data, gti=gti)
    return bak


def make_rsp(det_name):
    # 8 photon bins x 4 energy channels
    matrix = [[25.2, 0.0, 0.0, 0.0],
              [51.8, 54.9, 0.0, 0.0],
              [2.59, 82.0, 44.8, 0.0],
              [3.10, 11.6, 77.0, 0.13],
              [1.26, 6.21, 29.3, 14.6],
              [0.45, 3.46, 13.8, 9.98],
              [0.52, 4.39, 13.3, 3.93],
              [0.79, 7.14, 16.1, 3.92]]
    emin = [5.00, 15.8, 50.0, 158., 500., 1581, 5000, 15811]
    emax = [15.8, 50.0, 158., 500., 1581, 5000, 15811, 50000]
    chanlo = [4.60, 27.3, 102., 538.]
    chanhi = [27.3, 102., 538., 2000]
    drm = ResponseMatrix(matrix, emin, emax, chanlo, chanhi)

    tstart = 524666421.47
    tstop = 524666521.47
    trigtime = 524666471.47
    rsp = Rsp.from_data(drm, start_time=tstart, stop_time=tstop,
                        trigger_time=trigtime, detector=det_name)
    return rsp


class TestChisq(unittest.TestCase):

    def setUp(self):
        self.obs_counts = np.array([1, 10, 100, 0.0])
        self.back_rates = np.array([0.1, 1.0, 10.0, 0.0])
        self.back_var = np.array([0.01, 0.1, 1.0, 0.0])
        self.mod_rates = np.array([0.1, 1.0, 10.0, 0.0])
        self.exposure = 2.0

    def test_eval(self):
        f = chisq(self.obs_counts, self.back_rates, self.back_var,
                  self.mod_rates, self.exposure)
        self.assertAlmostEqual(f, -38.423, places=3)


class TestCstat(unittest.TestCase):

    def setUp(self):
        self.obs_counts = np.array([1, 10, 100, 0.0])
        self.back_counts = np.array([0.2, 2.0, 20.0, 0.0])
        self.mod_rates = np.array([0.1, 1.0, 10.0, 0.0])
        self.exposure = 2.0

    def test_eval(self):
        f = cstat(self.obs_counts, self.mod_rates, self.exposure,
                  self.back_counts, self.exposure)
        self.assertAlmostEqual(f, -18.709, places=3)


class TestPstat(unittest.TestCase):

    def setUp(self):
        self.obs_counts = np.array([1, 10, 100, 0.0])
        self.back_rates = np.array([0.1, 1.0, 10.0, 0.0])
        self.mod_rates = np.array([0.1, 1.0, 10.0, 0.0])
        self.exposure = 2.0

    def test_eval(self):
        f = pstat(self.obs_counts, self.mod_rates, self.exposure, self.back_rates)
        self.assertAlmostEqual(f, -35.108, places=3)


class TestPgstat(unittest.TestCase):

    def setUp(self):
        self.obs_counts = np.array([1, 10, 100, 0.0])
        self.back_counts = np.array([0.2, 2.0, 20.0, 0.0])
        self.back_var = np.array([0.01, 0.1, 1.0, 0.0])
        self.mod_rates = np.array([0.1, 1.0, 10.0, 0.0])
        self.exposure = 2.0

    def test_eval(self):
        f = pgstat(self.obs_counts, self.mod_rates, self.exposure,
                   self.back_counts, self.back_var, self.exposure)
        self.assertAlmostEqual(f, -33.931, places=3)


class TestSpectralFitterOne(unittest.TestCase):

    def setUp(self):
        self.pha = make_first_pha()
        self.bak = make_first_bak()
        self.rsp = make_rsp('det0')
        self.fitter = SpectralFitterChisq([self.pha], [self.bak.data],
                                          [self.rsp], method='SLSQP')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_covariance(self):
        covar = self.fitter.covariance.flatten()

        test_vals = [2.07e-6, 1.30e-5, 1.30e-5, -6.57e-4]
        for i in range(4):
            self.assertAlmostEqual(covar[i], test_vals[i], places=2)

    def test_detectors(self):
        self.assertListEqual(self.fitter.detectors, ['det0'])

    def test_dof(self):
        self.assertEqual(self.fitter.dof, 2)

    def test_energy_range(self):
        erange = self.fitter.energy_range
        self.assertTupleEqual(erange, (4.60, 2000.0))

    def test_function_components(self):
        self.assertIsNone(self.fitter.function_components)

    def test_function_name(self):
        self.assertEqual(self.fitter.function_name, 'PowerLaw')

    def test_hessian(self):
        hessian = self.fitter.hessian.flatten()
        test_vals = [-5.52657433e+5, -1.08974205e+4, -1.08974205e+4, -1.73589253e+3]
        npt.assert_allclose(hessian, test_vals)

    def test_jacobian(self):
        jac = self.fitter.jacobian.tolist()
        test_vals = [-3.29953594e-1, 7.17800293e-4]
        npt.assert_allclose(jac, test_vals)

    def test_message(self):
        self.assertEqual(self.fitter.message, 'Optimization terminated successfully')

    def test_num_components(self):
        self.assertEqual(self.fitter.num_components, 1)

    def test_num_sets(self):
        self.assertEqual(self.fitter.num_sets, 1)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_statistic(self):
        self.assertAlmostEqual(self.fitter.statistic, 2.98, places=2)

    def test_success(self):
        self.assertTrue(self.fitter.success)

    def test_symmetric_errors(self):
        errs = self.fitter.symmetric_errors
        test_vals = [0.0014, 0.0256]
        for i in range(2):
            self.assertAlmostEqual(errs[i], test_vals[i], places=4)

    def test_asymmetric_errors(self):
        errs = self.fitter.asymmetric_errors().flatten()
        test_vals = [0.002, 0.002, 0.036, 0.037]
        for i in range(4):
            self.assertAlmostEqual(errs[i], test_vals[i], places=3)

        errs = self.fitter.asymmetric_errors(cl=0.9).flatten()
        test_vals = [0.003, 0.003, 0.051, 0.052]
        for i in range(4):
            self.assertAlmostEqual(errs[i], test_vals[i], places=3)

        errs = self.fitter.asymmetric_errors(cl=0.99).flatten()
        test_vals = [0.004, 0.004, 0.071, 0.073]
        for i in range(4):
            self.assertAlmostEqual(errs[i], test_vals[i], places=3)

        with self.assertRaises(ValueError):
            self.fitter.asymmetric_errors(cl=-0.5)

        with self.assertRaises(ValueError):
            self.fitter.asymmetric_errors(cl=2.0)

    def test_data_count_spectrum(self):
        ecent, ewidths, counts, errs, ulmask = self.fitter.data_count_spectrum()
        self.assertListEqual(ecent[0].tolist(), self.pha.data.centroids.tolist())
        self.assertListEqual(ewidths[0][0, :].tolist(),
                             (self.pha.data.centroids - self.pha.data.lo_edges).tolist())
        self.assertListEqual(ewidths[0][1, :].tolist(),
                             (self.pha.data.hi_edges - self.pha.data.centroids).tolist())
        test_vals = [38.790, 13.713, 1.862, 0.042]
        for i in range(4):
            self.assertAlmostEqual(counts[0][i], test_vals[i], places=3)

        test_vals = [2.581, 0.899, 0.126, 0.014]
        for i in range(4):
            self.assertAlmostEqual(errs[0][i], test_vals[i], places=3)

    def test_model_count_spectrum(self):
        spec = self.fitter.model_count_spectrum()[0]
        self.assertListEqual(spec.lo_edges.tolist(), self.pha.data.lo_edges.tolist())
        self.assertListEqual(spec.lo_edges.tolist(), self.pha.data.lo_edges.tolist())
        self.assertListEqual(spec.exposure.tolist(), self.pha.data.exposure.tolist())

        test_vals = [215, 282, 194, 19]
        for i in range(4):
            self.assertAlmostEqual(spec.counts[i], test_vals[i], delta=1.0)

    def test_model_variance(self):
        var = self.fitter.model_variance()[0]
        test_vals = [6.660, 0.809, 0.016, 1.927e-4]
        for i in range(4):
            self.assertAlmostEqual(var[i], test_vals[i], places=3)

    def test_residuals(self):
        ecent, ewidths, resids, uncert = self.fitter.residuals()
        self.assertListEqual(ecent[0].tolist(), self.pha.data.centroids.tolist())
        self.assertListEqual(ewidths[0][0, :].tolist(),
                             (self.pha.data.centroids - self.pha.data.lo_edges).tolist())
        self.assertListEqual(ewidths[0][1, :].tolist(),
                             (self.pha.data.hi_edges - self.pha.data.centroids).tolist())

        test_vals = [6.87944898e-1, -1.12106660, 9.58694144e-1, -5.40612331e-1]
        npt.assert_allclose(resids[0], test_vals)

        self.assertListEqual(uncert[0].tolist(), [1.0] * 4)

        _, _, resids, uncerts = self.fitter.residuals(sigma=False)
        test_vals = [1.77541636, -1.00830332, 1.21027754e-1, -7.50400385e-3]
        npt.assert_allclose(resids[0], test_vals)

        test_vals = [2.58075373, 8.99414291e-1, 1.262423e-1, 1.38805636e-2]
        npt.assert_allclose(uncerts[0], test_vals)

    def test_sample_flux(self):
        fluxes = self.fitter.sample_flux((50.0, 300.0), num_samples=10)
        for i in range(10):
            self.assertAlmostEqual(fluxes[i], 8.5, delta=1.0)

    def test_sample_parameters(self):
        samples = self.fitter.sample_parameters(size=10)
        for i in range(10):
            self.assertAlmostEqual(samples[i, 0], 0.05, delta=0.03)
            self.assertAlmostEqual(samples[i, 1], -1.30, delta=0.3)

    def test_sample_spectrum(self):
        energies, func = self.fitter.sample_spectrum('photon', num_samples=10,
                                                     num_points=4)

        test_vals = [4.60, 34.85, 264.0, 2000.]
        for i in range(4):
            self.assertAlmostEqual(energies[i], test_vals[i], places=2)

        test_vals = [2.7, 2.0e-1, 1.4e-2, 1.0e-3]
        for i in range(10):
            for j in range(4):
                self.assertAlmostEqual(func[i, j] / test_vals[j], 1.0, delta=0.5)

        energies, func = self.fitter.sample_spectrum('energy', num_samples=10,
                                                     num_points=4)

        test_vals = [12.59, 6.86, 3.74, 2.04]
        for i in range(10):
            for j in range(4):
                self.assertAlmostEqual(func[i, j] / test_vals[j], 1.0, delta=0.5)

        energies, func = self.fitter.sample_spectrum('nufnu', num_samples=10,
                                                     num_points=4)

        test_vals = [57.93, 239.06, 986.49, 4070.91]
        for i in range(10):
            for j in range(4):
                self.assertAlmostEqual(func[i, j] / test_vals[j], 1.0, delta=0.5)

    def test_save_and_load(self):
        self.fitter.save('fitter.npz')
        fitter = SpectralFitterChisq.load('fitter.npz')
        os.remove('fitter.npz')
        self.assertListEqual(fitter.parameters.tolist(),
                             self.fitter.parameters.tolist())

    def test_spectrum(self):
        energies, func = self.fitter.spectrum('photon', num_points=4)

        test_vals = [4.60, 34.85, 264.0, 2000.]
        for i in range(4):
            self.assertAlmostEqual(energies[i], test_vals[i], places=2)

        test_vals = [2.7, 2.0e-1, 1.4e-2, 1.0e-3]
        for i in range(4):
            self.assertAlmostEqual(func[i] / test_vals[i], 1.0, delta=0.5)

        energies, func = self.fitter.spectrum('energy', num_points=4)

        test_vals = [12.59, 6.86, 3.74, 2.04]
        for i in range(4):
            self.assertAlmostEqual(func[i] / test_vals[i], 1.0, delta=0.5)

        energies, func = self.fitter.spectrum('nufnu', num_points=4)

        test_vals = [57.93, 239.06, 986.49, 4070.91]
        for i in range(4):
            self.assertAlmostEqual(func[i] / test_vals[i], 1.0, delta=0.5)

    def test_errors(self):

        with self.assertRaises(ValueError):
            self.fitter.asymmetric_errors(cl=-1.0)

        with self.assertRaises(ValueError):
            self.fitter.asymmetric_errors(cl=1.1)

        with self.assertRaises(ValueError):
            self.fitter.data_count_spectrum(upper_limits_sigma=-1.0)


class TestSpectralFitterTwo(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='TNC')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_detectors(self):
        self.assertListEqual(self.fitter.detectors, ['det0', 'det1'])

    def test_dof(self):
        self.assertEqual(self.fitter.dof, 6)

    def test_energy_range(self):
        erange = self.fitter.energy_range
        self.assertTupleEqual(erange, (4.60, 2000.0))

    def test_num_sets(self):
        self.assertEqual(self.fitter.num_sets, 2)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_statistic(self):
        self.assertAlmostEqual(self.fitter.statistic, 11.02, places=2)

    def test_success(self):
        self.assertTrue(self.fitter.success)

    def test_data_count_spectrum(self):
        ecent, ewidths, counts, errs, ulmask = self.fitter.data_count_spectrum()
        self.assertListEqual(ecent[0].tolist(), self.pha1.data.centroids.tolist())
        self.assertListEqual(ewidths[0][0, :].tolist(),
                             (self.pha1.data.centroids - self.pha1.data.lo_edges).tolist())
        self.assertListEqual(ewidths[0][1, :].tolist(),
                             (self.pha1.data.hi_edges - self.pha1.data.centroids).tolist())
        self.assertListEqual(ecent[1].tolist(), self.pha2.data.centroids.tolist())
        self.assertListEqual(ewidths[1][0, :].tolist(),
                             (self.pha2.data.centroids - self.pha2.data.lo_edges).tolist())
        self.assertListEqual(ewidths[1][1, :].tolist(),
                             (self.pha2.data.hi_edges - self.pha2.data.centroids).tolist())

        test_vals = [38.790, 13.713, 1.862, 0.042]
        for i in range(4):
            self.assertAlmostEqual(counts[0][i], test_vals[i], places=3)

        test_vals = [2.623, 0.910, 0.127, 0.014]
        for i in range(4):
            self.assertAlmostEqual(errs[0][i], test_vals[i], places=3)

        test_vals = [42.092, 14.060, 1.998, 0.036]
        for i in range(4):
            self.assertAlmostEqual(counts[1][i], test_vals[i], places=3)

        test_vals = [2.628, 0.911, 0.127, 0.014]
        for i in range(4):
            self.assertAlmostEqual(errs[1][i], test_vals[i], places=3)

    def test_model_count_spectrum(self):
        spec = self.fitter.model_count_spectrum()
        self.assertListEqual(spec[0].lo_edges.tolist(), self.pha1.data.lo_edges.tolist())
        self.assertListEqual(spec[0].lo_edges.tolist(), self.pha1.data.lo_edges.tolist())
        self.assertListEqual(spec[0].exposure.tolist(), self.pha1.data.exposure.tolist())

        test_vals = [222, 289, 197, 19]
        for i in range(4):
            self.assertAlmostEqual(spec[0].counts[i], test_vals[i], delta=1.0)

        self.assertListEqual(spec[1].lo_edges.tolist(), self.pha2.data.lo_edges.tolist())
        self.assertListEqual(spec[1].lo_edges.tolist(), self.pha2.data.lo_edges.tolist())
        self.assertListEqual(spec[1].exposure.tolist(), self.pha2.data.exposure.tolist())

        test_vals = [222, 289, 197, 19]
        for i in range(4):
            self.assertAlmostEqual(spec[1].counts[i], test_vals[i], delta=1.0)

    def test_model_variance(self):
        var = self.fitter.model_variance()

        test_vals = [6.881, 0.829, 0.016, 1.927e-4]
        for i in range(4):
            self.assertAlmostEqual(var[0][i], test_vals[i], places=3)

        test_vals = [6.905, 0.830, 0.016, 1.969e-4]
        for i in range(4):
            self.assertAlmostEqual(var[1][i], test_vals[i], places=3)

    def test_residuals(self):
        ecent, ewidths, resids, uncert = self.fitter.residuals()
        self.assertListEqual(ecent[0].tolist(), self.pha1.data.centroids.tolist())
        self.assertListEqual(ewidths[0][0, :].tolist(),
                             (self.pha1.data.centroids - self.pha1.data.lo_edges).tolist())
        self.assertListEqual(ewidths[0][1, :].tolist(),
                             (self.pha1.data.hi_edges - self.pha1.data.centroids).tolist())
        test_vals = [0.187, -1.519, 0.763, -0.539]
        for i in range(4):
            self.assertAlmostEqual(resids[0][i], test_vals[i], places=3)
        self.assertListEqual(uncert[0].tolist(), [1.0] * 4)

        self.assertListEqual(ecent[1].tolist(), self.pha2.data.centroids.tolist())
        self.assertListEqual(ewidths[1][0, :].tolist(),
                             (self.pha2.data.centroids - self.pha2.data.lo_edges).tolist())
        self.assertListEqual(ewidths[1][1, :].tolist(),
                             (self.pha2.data.hi_edges - self.pha2.data.centroids).tolist())
        test_vals = [1.444, -1.138, 1.834, -1.028]
        for i in range(4):
            self.assertAlmostEqual(resids[1][i], test_vals[i], places=3)
        self.assertListEqual(uncert[1].tolist(), [1.0] * 4)

        _, _, resids, uncerts = self.fitter.residuals(sigma=False)
        test_vals = [0.491, -1.382, 0.097, -0.007]
        for i in range(4):
            self.assertAlmostEqual(resids[0][i], test_vals[i], places=3)
        test_vals = [2.623, 0.910, 0.127, 0.014]
        for i in range(4):
            self.assertAlmostEqual(uncerts[0][i], test_vals[i], places=3)

        test_vals = [3.793, -1.036, 0.234, -0.014]
        for i in range(4):
            self.assertAlmostEqual(resids[1][i], test_vals[i], places=3)
        test_vals = [2.628, 0.911, 0.127, 0.014]
        for i in range(4):
            self.assertAlmostEqual(uncerts[1][i], test_vals[i], places=3)


class TestNelderMead(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='Nelder-Mead')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestPowell(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='Powell')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestConjugateGradient(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='CG')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl, options={'gtol': 0.5})

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestBFGS(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='BFGS')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl, options={'gtol': 0.1})

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestNewton(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='Newton-CG')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestLBFGSB(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='L-BFGS-B')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestCOBYLA(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='COBYLA')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestTrustConstr(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='trust-constr')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestTrustNCG(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='trust-ncg')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestTrustKrylov(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='trust-krylov')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestDogleg(unittest.TestCase):
    def test_fit(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='dogleg')
        pl = PowerLaw()
        pl.max_values[1] = 10.0

        with self.assertRaises(RuntimeError):
            self.fitter.fit(pl)


class TestTrustExact(unittest.TestCase):
    def test_fit(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='trust-exact')
        pl = PowerLaw()
        pl.max_values[1] = 10.0

        with self.assertRaises(RuntimeError):
            self.fitter.fit(pl)


class TestSpectralFitterCstat(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterCstat([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='TNC')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestSpectralFitterPstat(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterPstat([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='TNC')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestSpectralFitterPgstat(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterPgstat([self.pha1, self.pha2],
                                           [self.bak1.data, self.bak2.data],
                                           [self.rsp1, self.rsp2], method='TNC')
        pl = PowerLaw()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl)

    def test_parameters(self):
        self.assertAlmostEqual(self.fitter.parameters[0], 0.05, places=2)
        self.assertAlmostEqual(self.fitter.parameters[1], -1.3, places=1)

    def test_success(self):
        self.assertTrue(self.fitter.success)


class TestFitTwoComponents(unittest.TestCase):
    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp1 = make_rsp('det0')
        self.rsp2 = make_rsp('det1')
        self.fitter = SpectralFitterChisq([self.pha1, self.pha2],
                                          [self.bak1.data, self.bak2.data],
                                          [self.rsp1, self.rsp2], method='TNC')
        pl = PowerLaw()
        bb = BlackBody()
        pl.max_values[1] = 10.0
        self.fitter.fit(pl + bb)

    def test_success(self):
        self.assertTrue(self.fitter.success)

    def test_covariance(self):
        covar = self.fitter.covariance
        self.assertTupleEqual(covar.shape, (4, 4))

    def test_dof(self):
        self.assertEqual(self.fitter.dof, 4)

    def test_function_components(self):
        self.assertListEqual(self.fitter.function_components,
                             ['PowerLaw', 'BlackBody'])

    def test_function_name(self):
        self.assertEqual(self.fitter.function_name, 'PowerLaw + BlackBody')

    def test_num_components(self):
        self.assertEqual(self.fitter.num_components, 2)

    def test_sample_spectrum(self):
        spec = self.fitter.sample_spectrum('counts', num_samples=1,
                                           num_points=4, components=True)
        self.assertEqual(len(spec), 2)

    def test_spectrum(self):
        spec = self.fitter.spectrum('counts', num_points=4, components=True)
        self.assertEqual(len(spec), 2)


class TestFailures(unittest.TestCase):

    def setUp(self):
        self.pha1 = make_first_pha()
        self.bak1 = make_first_bak()
        self.pha2 = make_second_pha()
        self.bak2 = make_second_bak()
        self.rsp = make_rsp('det0')
        self.fitter = SpectralFitterChisq([self.pha1], [self.bak1.data], [self.rsp])

    def test_wrong_num_inputs(self):
        with self.assertRaises(ValueError):
            SpectralFitterChisq([self.pha1, self.pha2], [self.bak1.data],
                                [self.rsp, self.rsp])

    def test_wrong_length_channel_masks(self):
        chan_mask = np.ones(4, dtype=bool)
        with self.assertRaises(ValueError):
            SpectralFitterChisq([self.pha1], [self.bak1.data], [self.rsp],
                                channel_masks=[chan_mask, chan_mask])

    def test_wrong_background(self):
        with self.assertRaises(ValueError):
            SpectralFitterChisq([self.pha1], [self.bak1.data.counts], [self.rsp])

    def test_incorrect_asymmetric_errors_input(self):
        with self.assertRaises(RuntimeError):
            self.fitter.asymmetric_errors()

    def test_model_count_spectrum(self):
        with self.assertRaises(RuntimeError):
            self.fitter.model_count_spectrum()

    def test_model_variance(self):
        with self.assertRaises(RuntimeError):
            self.fitter.model_variance()

    def test_residuals(self):
        with self.assertRaises(RuntimeError):
            self.fitter.residuals()

    def test_sample_flux(self):
        with self.assertRaises(RuntimeError):
            self.fitter.sample_flux((50.0, 300.))

    def test_sample_parameters(self):
        with self.assertRaises(RuntimeError):
            self.fitter.sample_parameters()

    def test_sample_spectrum(self):
        with self.assertRaises(RuntimeError):
            self.fitter.sample_spectrum('counts')

    def test_spectrum(self):
        with self.assertRaises(RuntimeError):
            self.fitter.spectrum('counts')
