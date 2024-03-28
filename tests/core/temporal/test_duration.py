#
#     Authors: Oliver J. Roberts (USRA)
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
from gdt.core.temporal.duration import Duration
import numpy as np
from gdt.core.background.primitives import BackgroundSpectrum
from gdt.core.simulate.profiles import linear
from gdt.core.simulate.tte import TteBackgroundSimulator
from gdt.core.spectra.functions import Band
from gdt.core.simulate.profiles import norris
from gdt.core.simulate.tte import TteSourceSimulator
from gdt.core.data_primitives import ResponseMatrix
from gdt.core.response import Rsp
from gdt.core.tte import PhotonList
from gdt.core.binning.unbinned import bin_by_time
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial

rates = [17.5, 77.5, 37.5, 57.5]
rate_uncert = [1.896, 1.889, 1.919, 1.66]
emin = [4.60, 27.3, 102., 538.]
emax = [27.3, 102., 538., 2000]
exposure = 0.128
back_spec = BackgroundSpectrum(rates, rate_uncert, emin, emax, exposure)
back_sim = TteBackgroundSimulator(back_spec, 'Gaussian', linear, (40.0, 0.1), deadtime=1e-6)

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
rsp = Rsp.from_data(drm, start_time=tstart, stop_time=tstop, trigger_time=trigtime, detector='det0')

# (amplitude, Epeak, alpha, beta)
band_params = (0.1, 567.0, -0.7, -3.2)
# (amplitude, tstart, trise, tdecay)
norris_params = (1.5, 1.47, 0.5, 1.0)
norris_params2 = (0.5, 3.47, 0.9, 0.3)
norris_params3 = (0.8, -0.47, 0.5, 0.2)
src_sim = TteSourceSimulator(rsp, Band(), band_params, norris, norris_params,
                             deadtime=1e-6)
src_sim2 = TteSourceSimulator(rsp, Band(), band_params, norris, norris_params2,
                              deadtime=1e-6)
src_sim3 = TteSourceSimulator(rsp, Band(), band_params, norris, norris_params3,
                              deadtime=1e-6)

back_tte = back_sim.to_tte(-10.0, 30.0)
src_tte = src_sim.to_tte(-10.0, 30.0)
src_tte2 = src_sim2.to_tte(-10.0, 30.0)
src_tte3 = src_sim3.to_tte(-10.0, 30.0)
total_tte = PhotonList.merge([back_tte, src_tte, src_tte2, src_tte3])

# src FILES INPUT ################
duration_interval = (0.05, 0.95)
num_sims = 10000
confidence = 0.9
time_Res = 0.016  # in s
nai_erange = (10.0, 1000.0)
nai_50_300 = (50.0, 300.0)
view_range = (-10, 20)  # zoom in to this time range
bkgd_range = [(-10, -2), (10, 20)]  # the background fit ranges

tte = total_tte.to_phaii(bin_by_time, time_Res, time_ref=0.0).slice_time(view_range).slice_energy(nai_erange)
phaii = tte.data.integrate_energy(nai_erange[0], nai_erange[1])
timebins_list = [phaii]

# Background FILES INPUT ################

bf_phaii = BackgroundFitter.from_phaii(tte, Polynomial, time_ranges=bkgd_range)
bf_phaii.fit(order=2)
bkgds_phaiis = bf_phaii.interpolate_bins(tte.data.tstart, tte.data.tstop).integrate_energy(nai_erange[0],
                                                                                           nai_erange[1])
bkgds_list = [bkgds_phaiis]

timebins_x = timebins_list[0].centroids
timebins_y = timebins_list[0].counts
br = bkgds_list[0].counts

workingrate = timebins_y.T[:, ] - br.T[:, ]
cumflsum = workingrate.cumsum(axis=1)

array = (timebins_x, cumflsum.T[:, ])
dur_per = duration_interval[1] - duration_interval[0]
tbins = timebins_x



class TestDuration(unittest.TestCase):

    def setUp(self):
        self.duration = Duration(timebins_list, bkgds_list, duration_interval)

    def test_inputExists(self):
        self.assertIsNotNone(self.duration)

    def test_inputListEq(self):
        self.assertEqual(len(self.duration.timebins_list), len(self.duration.bkgds_list))

    def test_inputType(self):
        self.assertIsInstance(self.duration.duration_interval[0], float)
        self.assertIsInstance(self.duration.duration_interval[1], float)
        self.assertIsInstance(self.duration.error_prop(2, 2), float)
        self.assertEqual(len(self.duration.timebins_list[0].counts), len(self.duration.bkgds_list[0].counts))
        self.assertAlmostEqual(self.duration.timebins_list[0].centroids[0], self.duration.bkgds_list[0].time_range[0], places=1)
        self.assertAlmostEqual(self.duration.timebins_list[0].centroids[-1], self.duration.bkgds_list[0].time_range[-1], places=1)
        self.assertIsInstance(duration_interval[0], float)
        self.assertIsInstance(duration_interval[1], float)

    def test_functionReturn(self):
        self.assertLessEqual(self.duration.duration_interval[0], 1.0)
        self.assertLessEqual(self.duration.duration_interval[1], 1.0)
        self.assertIsInstance(self.duration.calculate(num_sims, confidence)[0], float)
        self.assertIsInstance(self.duration.calculate(num_sims, confidence)[1], float)
        self.assertIsInstance(self.duration.calculate(num_sims, confidence)[2], float)
        self.assertIsNotNone(self.duration.calculate(num_sims, confidence)[0])
        self.assertIsNotNone(self.duration.calculate(num_sims, confidence)[1])
        self.assertIsNotNone(self.duration.calculate(num_sims, confidence)[2])
        self.assertEqual(self.duration.error_prop(2, 2), np.sqrt(8))
        self.assertIsInstance(self.duration.findtparams(array[1], dur_per, tbins), float)
        self.assertIsInstance(self.duration.quantiles(self.duration.findtparams(array[1], dur_per, tbins),0.9)[0], float)
        self.assertIsInstance(self.duration.quantiles(self.duration.findtparams(array[1], dur_per, tbins), 0.9)[1], float)

    def test_inputListInp(self):
        self.assertEqual(len(self.duration.timebins_list[0].centroids), len(self.duration.bkgds_list[0].time_centroids))
        self.assertEqual(len(self.duration.timebins_list[0].range), len(self.duration.bkgds_list[0].time_range))

if __name__ == '__main__':
    unittest.main()