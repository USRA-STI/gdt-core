#
#     Authors: Oliver Roberts (USRA)
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
import numpy as np
import os
import unittest
from unittest import TestCase
from gdt.core.duration import Duration
import gdt.core
import numpy as np
from gdt.core import data_path
from gdt.missions.fermi.gbm.collection import GbmDetectorCollection
from gdt.missions.fermi.gbm.tte import GbmTte
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial
from gdt.core.binning.unbinned import bin_by_time
from gdt.core.data_primitives import TimeBins


filepath = data_path.joinpath('/Users/orobert2/PycharmProjects/gdt-core_roberts_dev/tests/core/test_data/')

#inputs:
bn = 'bn230812790'

n0 = GbmTte.open(str(filepath) + '/' + str(bn) + '/glg_tte_n0_{}_v04.fit'.format(bn))
n6 = GbmTte.open(str(filepath) + '/' + str(bn) + '/glg_tte_n6_{}_v04.fit'.format(bn))
n7 = GbmTte.open(str(filepath) + '/' + str(bn) + '/glg_tte_n7_{}_v04.fit'.format(bn))

time_Res = 0.016 #in s
nai_erange = (10.0, 1000.0)
nai_50_300 = (50.0, 300.0)
view_range = (-10, 25) # zoom in to this time range
bkgd_range = [(-10, -0.5),(10, 25)] # the background fit ranges

tte_n0 = n0.to_phaii(bin_by_time, time_Res, time_ref=0.0).slice_time(view_range).slice_energy(nai_erange)
tte_n6 = n6.to_phaii(bin_by_time, time_Res, time_ref=0.0).slice_time(view_range).slice_energy(nai_erange)
tte_n7 = n7.to_phaii(bin_by_time, time_Res, time_ref=0.0).slice_time(view_range).slice_energy(nai_erange)

phaii_n0 = tte_n0.data.integrate_energy(nai_erange[0],nai_erange[1])
phaii_n6 = tte_n6.data.integrate_energy(nai_erange[0],nai_erange[1])
phaii_n7 = tte_n7.data.integrate_energy(nai_erange[0],nai_erange[1])

phaiis = GbmDetectorCollection.from_list([tte_n0,tte_n6,tte_n7], names=['n0','n6','n7'], dets=['n0','n6','n7'])

# timebins = [phaii_n0,phaii_n6,phaii_n7]
# timebins = [phaii_n0]


################ Background ################

bf_phaiis = [BackgroundFitter.from_phaii(phaii, Polynomial, time_ranges=bkgd_range) for phaii in phaiis]
for bf_phaii in bf_phaiis:
    bf_phaii.fit(order=2)

bk_list=[]
for bf_phaii in bf_phaiis:
    bkgds_phaiis = bf_phaii.interpolate_bins(phaiis.data()[0].tstart, phaiis.data()[0].tstop).integrate_energy(nai_erange[0],nai_erange[1])
    bk_list.append(bkgds_phaiis)


duration_interval = (0.25,0.75)
num_sims = 10000
confidence = 0.67

# timebins_list = TimeBins.sum([phaii_n0,phaii_n6,phaii_n7])
# data = TimeBins.sum([phaii_n0, phaii_n6, phaii_n7])
# timebins_list_x = timebins_list.centroids
# timebins_list_y = timebins_list.counts


class TestDuration(unittest.TestCase):

    timebins = TimeBins([phaii_n0,phaii_n6,phaii_n7])
    duration = Duration(timebins, bk_list, duration_interval)
    duration.calculate(num_sims,confidence)
    data = TimeBins.sum([phaii_n0, phaii_n6, phaii_n7])

    def test_inputExists(self):
        self.assertIsNotNone(self.duration)

    # def test_data(self):
    #     self.asserttrue(self.duration.timebins_list, TimeBins)
    #
    # # def test_attributes(self):
    # #     self.assertEqual(len(self.duration), 3)


        
if __name__ == '__main__':
    unittest.main()