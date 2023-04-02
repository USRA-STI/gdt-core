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
import numpy.testing as npt
from gdt.core.spectra.functions import *


class MyFirstFunction(Function):
    nparams = 2
    param_list = [('C0', 'units1', 'Coeff1'),
                  ('C1', 'units2', 'Coeff2')]
    default_values = [0.0, 0.0]
    delta_abs = [0.1, 0.1]
    delta_rel = [0.01, 0.01]
    min_values = [-np.inf, -np.inf]
    max_values = [np.inf, np.inf]
    free = [True, True]

    def eval(self, params, x):
        return params[0] + params[1] * x


class MySecondFunction(Function):
    nparams = 2
    param_list = [('C0', 'units1', 'Coeff1'),
                  ('C1', 'units2', 'Coeff2')]
    default_values = [0.0, 0.0]
    delta_abs = [0.1, 0.1]
    delta_rel = [0.01, 0.01]
    min_values = [-np.inf, -np.inf]
    max_values = [np.inf, np.inf]
    free = [True, True]

    def eval(self, params, x):
        return params[0] + params[1] * x ** 2


class MyThirdFunction(Function):
    nparams = 1

    def eval(self, params, x):
        arr = np.empty(x.size)
        arr.fill(params[0])
        return arr


class MyIncompleteFunction(Function):
    def eval(self, params, x):
        pass


class TestFunction(unittest.TestCase):

    def test_attributes(self):
        myfunc = MyFirstFunction()
        self.assertEqual(myfunc.name, 'MyFirstFunction')
        self.assertEqual(myfunc.nparams, 2)
        self.assertEqual(myfunc.num_components, 1)
        self.assertEqual(myfunc.param_list,
                         [('C0', 'units1', 'Coeff1'), ('C1', 'units2', 'Coeff2')])
        self.assertEqual(myfunc.free, [True, True])
        myfunc.free = [False, True]
        self.assertEqual(myfunc.free, [False, True])

    def test_eval(self):
        myfunc = MyFirstFunction()
        x = np.array([0.0, 1.0, 2.0])
        params = [10.0, 1.0]
        y = myfunc.eval(params, x)
        self.assertCountEqual(y, np.array([10.0, 11.0, 12.0]))

    def test_fit_eval(self):
        myfunc = MyFirstFunction()
        myfunc.free = [False, True]
        myfunc.default_values = [-10.0, 0.0]

        params = [1.0]
        x = np.array([0.0, 1.0, 2.0])
        y = myfunc.fit_eval(params, x)
        self.assertCountEqual(y, np.array([-10.0, -9.0, -8.0]))

    def test_parameter_bounds(self):
        myfunc = MyFirstFunction()
        myfunc.free = [False, True]
        bounds = myfunc.parameter_bounds(apply_state=False)
        self.assertEqual(bounds, [(-np.inf, np.inf), (-np.inf, np.inf)])
        bounds = myfunc.parameter_bounds(apply_state=True)
        self.assertEqual(bounds, [(-np.inf, np.inf)])

    def test_integrate(self):
        myfunc = MyFirstFunction()
        params = [10.0, 1.0]
        integral = myfunc.integrate(params, (1.0, 10.0), log=False, energy=False)
        self.assertEqual(integral, 139.5)
        integral = myfunc.integrate(params, (1.0, 10.0), log=True, energy=False)
        self.assertEqual(integral, 139.5)
        integral = myfunc.integrate(params, (1.0, 10.0), log=False, energy=True)
        self.assertAlmostEqual(integral, 828.0 * 1.602e-9)

        myfunc.free = [False, True]
        integral = myfunc.integrate([1.0], (1.0, 10.0), log=False, energy=False)
        self.assertAlmostEqual(integral, 49.5)

    def test_repr(self):
        myfunc = MyFirstFunction()
        expected_repr = ('<MyFirstFunction: 2 parameters;\n'
                         ' Defaults: C0 = 0.0 units1;\n'
                         '           C1 = 0.0 units2>')
        self.assertEqual(repr(myfunc), expected_repr)

    def test_incomplete(self):
        with self.assertRaises(AttributeError):
            MyIncompleteFunction()

    def test_incorrect_number_parameters(self):
        myfunc = MyFirstFunction()
        params = [10.0, 1.0, 100.0, 1000.0]
        with self.assertRaises(ValueError):
            myfunc.integrate(params, (1.0, 10.0), log=False, energy=False)


class TestAdditiveSuperFunction(unittest.TestCase):

    def test_attributes(self):
        myfunc = MyFirstFunction() + MySecondFunction()
        self.assertEqual(myfunc.name, 'MyFirstFunction + MySecondFunction')
        self.assertEqual(myfunc.nparams, 4)
        self.assertEqual(myfunc.num_components, 2)
        self.assertEqual(myfunc.param_list,
                         [('MyFirstFunction: C0', 'units1', 'Coeff1'),
                          ('MyFirstFunction: C1', 'units2', 'Coeff2'),
                          ('MySecondFunction: C0', 'units1', 'Coeff1'),
                          ('MySecondFunction: C1', 'units2', 'Coeff2')])
        self.assertEqual(myfunc.free, [True] * 4)

    def test_eval(self):
        myfunc = MyFirstFunction() + MySecondFunction()
        x = np.array([0.0, 1.0, 2.0])
        params = [10.0, 1.0, 10.0, 2.0]
        y = myfunc.eval(params, x)
        self.assertCountEqual(y, np.array([20.0, 23.0, 30.0]))

        y = myfunc.eval(params, x, components=True)
        self.assertCountEqual(y[0], np.array([10.0, 11.0, 12.0]))
        self.assertCountEqual(y[1], np.array([10.0, 12.0, 18.0]))


class TestMultiplicativeSuperFunction(unittest.TestCase):

    def test_attributes(self):
        myfunc = MyFirstFunction() * MySecondFunction()
        self.assertEqual(myfunc.name, 'MyFirstFunction * MySecondFunction')
        self.assertEqual(myfunc.nparams, 4)
        self.assertEqual(myfunc.num_components, 2)
        self.assertEqual(myfunc.param_list,
                         [('MyFirstFunction: C0', 'units1', 'Coeff1'),
                          ('MyFirstFunction: C1', 'units2', 'Coeff2'),
                          ('MySecondFunction: C0', 'units1', 'Coeff1'),
                          ('MySecondFunction: C1', 'units2', 'Coeff2')])
        self.assertEqual(myfunc.free, [True] * 4)

    def test_eval(self):
        myfunc = MyFirstFunction() * MySecondFunction()
        x = np.array([0.0, 1.0, 2.0])
        params = [10.0, 1.0, 10.0, 2.0]
        y = myfunc.eval(params, x)
        self.assertCountEqual(y, np.array([100., 132., 216.]))

        y = myfunc.eval(params, x, components=True)
        self.assertCountEqual(y[0], np.array([10.0, 11.0, 12.0]))
        self.assertCountEqual(y[1], np.array([10.0, 12.0, 18.0]))


class TestMixedMultipleSuperFunctions(unittest.TestCase):
    def test_attributes(self):
        myfunc = (MyFirstFunction() + MySecondFunction()) * MyThirdFunction()
        self.assertEqual(myfunc.name, 'MyFirstFunction + MySecondFunction * MyThirdFunction')
        self.assertEqual(myfunc.nparams, 5)
        self.assertEqual(myfunc.num_components, 3)

    def test_eval(self):
        myfunc = MyFirstFunction() + MySecondFunction() * MyThirdFunction()
        x = np.array([0.0, 1.0, 2.0])
        params = [10.0, 1.0, 10.0, 2.0, 100.0]
        y = myfunc.eval(params, x)
        self.assertCountEqual(y, np.array([2000.0, 2300.0, 3000.0]))

        y = myfunc.eval(params, x, components=True)
        self.assertCountEqual(y[0], np.array([10.0, 11.0, 12.0]))
        self.assertCountEqual(y[1], np.array([10.0, 12.0, 18.0]))
        self.assertCountEqual(y[2], np.array([100.0, 100.0, 100.0]))


class TestSpectralFunctions(unittest.TestCase):
    energies = np.array([10.0, 100.0, 1000.0])

    def test_powerlaw(self):
        params = (0.1, -2.0, 100.0)
        y = PowerLaw().eval(params, self.energies)
        npt.assert_array_equal(y, np.array([10.0, 0.1, 0.001]))

        params = (0.3, -1.5, 320)
        y = PowerLaw().eval(params, self.energies)
        expected = np.array([5.43058008E+01, 1.71730021E+00, 5.43058008E-02])
        npt.assert_allclose(y, expected)

    def test_comptonized(self):
        params = (0.1, 200.0, -1.0, 100.0)
        y = Comptonized().eval(params, self.energies)
        expected = np.array((9.51229425E-01, 6.06530660E-02, 6.73794700E-05))
        npt.assert_allclose(y, expected)

    def test_comptonized_bad_index(self):
        params = (0.1, 200.0, -2.0, 100.0)
        y = Comptonized().eval(params, self.energies)
        expected = np.zeros_like(self.energies)
        npt.assert_array_equal(expected, y)

    def test_band(self):
        params = (0.1, 200.0, -1.0, -2.3, 100.0)
        y = Band().eval(params, self.energies)
        true = np.array((0.951, 0.061, 4.73e-4))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_band_old(self):
        params = (0.1, 200.0, -1.0, -2.3, 100.0)
        y = BandOld().eval(params, self.energies)
        true = np.array((0.951, 0.061, 4.73e-4))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_bpl(self):
        params = (0.01, 700.0, -1.0, -2.0, 100.0)
        y = BrokenPowerLaw().eval(params, self.energies)
        true = np.array((0.100, 0.010, 7e-4))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_2bpl(self):
        params = (0.01, 90.0, 500, -0.5, -1.0, -2.0, 100.0)
        y = DoubleBrokenPowerLaw().eval(params, self.energies)
        true = np.array((0.032, 0.009, 0.015))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_smoothly_bpl(self):
        params = (0.01, 700.0, 0.3, -1.0, -2.0, 100.0)
        y = SmoothlyBrokenPowerLaw().eval(params, self.energies)
        true = np.array((0.100, 0.010, 6.3e-4))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_lognormal(self):
        params = (0.1, 5.0, 1.0)
        y = LogNormal().eval(params, self.energies)
        true = np.array((1.0e-4, 3.7e-4, 6.5e-6))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_gauss_log(self):
        params = (0.1, 100.0, 1.0)
        y = GaussianLog().eval(params, self.energies)
        true = np.array((0.006, 0.094, 0.006))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_gauss_log_fwhm(self):
        params = (0.1, 100.0, 1.0, 0.1)
        y = GaussianLogVaryingFWHM().eval(params, self.energies)
        true = np.array((0.003, 0.094, 0.009))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_sunyaev_titarchuk(self):
        params = (0.1, 30.0, 10.0, 3.0)
        y = SunyaevTitarchuk().eval(params, self.energies)
        true = np.array((0.003, 2e-40, 0.00))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_ottb(self):
        params = (0.1, 30.0, 100.0)
        y = OTTB().eval(params, self.energies)
        true = np.array((20.086, 0.100, 9.4e-16))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_blackbody(self):
        params = (0.1, 30.0)
        y = BlackBody().eval(params, self.energies)
        true = np.array((25.277, 36.994, 3.3e-10))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_yang_soong(self):
        params = (0.1, -1.0, 200.0, 300.0, 50.0, 200.0, 50.0, 100.0)
        y = YangSoongPulsar().eval(params, self.energies)
        true = np.array((0.010, 0.100, 3.5e-6))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_tanaka(self):
        params = (0.10, -1.0, 200.0, 1.0, 1.0, 1, 100.0, 50.0, 100.0)
        y = TanakaPulsar().eval(params, self.energies)
        true = np.array((0.004, 0.010, 0.002))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_tanaka_NL0(self):
        params = (0.10, -1.0, 200.0, 1.0, 1.0, 0, 100.0, 50.0, 100.0)
        y = TanakaPulsar().eval(params, self.energies)
        test_vals = [1.56831219e-3, 0.01, 1.11089965e-03]
        npt.assert_allclose(y, test_vals)

    def test_otts(self):
        params = (0.1, 100.0)
        y = OTTS().eval(params, self.energies)
        true = np.array((0.159, 0.272, 0.862))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_gaussline(self):
        params = (0.1, 100.0, 8.0)
        y = GaussLine().eval(params, self.energies)
        true = np.array((5.6e-4, 5.6e-5, 0.0))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_low_cutoff(self):
        params = (100.0, 10.0)
        y = LowEnergyCutoff().eval(params, self.energies)
        true = np.array((8.1e-7, 1.0, 1.0))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_high_cutoff(self):
        params = (100.0, 100.0)
        y = HighEnergyCutoff().eval(params, self.energies)
        true = np.array((1., 1., 0.001))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_multiplicative_pl(self):
        params = (-2.0, 100.0)
        y = PowerLawMult().eval(params, self.energies)
        true = np.array((100., 1.0, 0.01))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_multiplicative_gaussline(self):
        params = (1.0, 10.0, 8.0)
        y = GaussLineMult().eval(params, self.energies)
        true = np.array((1.125, 1.0, 1.0))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]

    def test_multiplicative_lorentzline(self):
        params = (1.0, 10.0, 8.0)
        y = LorentzLineMult().eval(params, self.energies)
        true = np.array((1.016, 1.000, 1.0))
        [self.assertAlmostEqual(y[i], true[i], places=3) for i in range(3)]
