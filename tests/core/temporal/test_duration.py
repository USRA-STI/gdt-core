# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 425 3rd Street SW, Suite 950, Washington, DC 20024
#
# Copyright 2017-2024 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland, Adam Goldstein and Oliver Roberts
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
from pathlib import Path
import numpy as np
import os
import unittest
from unittest import TestCase
from gdt.core.temporal.duration import Duration
import gdt.core
import numpy as np
from gdt.missions.fermi.gbm.collection import GbmDetectorCollection
from gdt.missions.fermi.gbm.tte import GbmTte
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial
from gdt.core.binning.unbinned import bin_by_time
from gdt.core.data_primitives import TimeBins

data_path = Path(__file__).parent.joinpath('data')
n0_file = data_path / 'glg_tte_n0_bn230812790_v04.fit'
n6_file = data_path / 'glg_tte_n6_bn230812790_v04.fit'
n7_file = data_path / 'glg_tte_n7_bn230812790_v04.fit'

n0 = GbmTte.open(n0_file)
n6 = GbmTte.open(n6_file)
n7 = GbmTte.open(n7_file)

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


duration_interval = (0.25, 0.75)
num_sims = 10000
confidence = 0.67


class TestDuration(unittest.TestCase):
    
    timebins = [phaii_n0, phaii_n6, phaii_n7]
    duration = Duration(timebins, bk_list, duration_interval)
    
    output = duration.calculate(num_sims, confidence)
    data = TimeBins.sum([phaii_n0, phaii_n6, phaii_n7])

    def test_inputExists(self):
        self.assertIsNotNone(self.output)

