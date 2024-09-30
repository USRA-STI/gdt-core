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
import unittest

from gdt.core.fits import exponential_card, fixed_card

class TestFitCards(unittest.TestCase):
    def test_fixed_card(self):
        val = 7.428703703703703e-4
        card = fixed_card('FIXEDV', val, comment='fixed floating point')
        assert card.image == "FIXEDV  = 0.00074              / fixed floating point                           "

    def test_exponential_card_float(self):
        val = 7.428703703703703e-4
        card = exponential_card('FLOATV', val, comment='exponential floating point')
        assert card.image == "FLOATV  = 7.42870E-4           / exponential floating point                     "

    def test_exponential_card_double(self):
        val = 7.428703703703703e-4
        card = exponential_card('DOUBLV', val, places=15, use_double=True,
                                comment='exponential floating point (as double)')
        assert card.image == "DOUBLV  = 7.428703703703703D-4 / exponential floating point (as double)         "
