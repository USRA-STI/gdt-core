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
import numpy as np
from gdt.core.detector import Detectors

class MyDetector(Detectors):
    det0 = ('Det0', 0,  0.0, 15.0)
    det1 = ('Det1', 1, 15.0, 30.0)
    det2 = ('Det2', 2, 30.0, 45.0)
    det3 = ('Det3', 3, 45.0, 60.0)

class TestDetector(unittest.TestCase):
    
    def test_number(self):
        self.assertEqual(MyDetector.det0.number, 0)
        self.assertEqual(MyDetector.det1.number, 1)
        self.assertEqual(MyDetector.det2.number, 2)
        self.assertEqual(MyDetector.det3.number, 3)

    def test_name(self):
        self.assertEqual(MyDetector.det0.name, 'det0')
        self.assertEqual(MyDetector.det1.name, 'det1')
        self.assertEqual(MyDetector.det2.name, 'det2')
        self.assertEqual(MyDetector.det3.name, 'det3')

    def test_full_name(self):
        self.assertEqual(MyDetector.det0.full_name, 'Det0')
        self.assertEqual(MyDetector.det1.full_name, 'Det1')
        self.assertEqual(MyDetector.det2.full_name, 'Det2')
        self.assertEqual(MyDetector.det3.full_name, 'Det3')

    def test_azimuth(self):
        self.assertEqual(MyDetector.det0.azimuth, 0.0)
        self.assertEqual(MyDetector.det1.azimuth, 15.0)
        self.assertEqual(MyDetector.det2.azimuth, 30.0)
        self.assertEqual(MyDetector.det3.azimuth, 45.0)

    def test_zenith(self):
        self.assertEqual(MyDetector.det0.zenith, 15.0)
        self.assertEqual(MyDetector.det1.zenith, 30.0)
        self.assertEqual(MyDetector.det2.zenith, 45.0)
        self.assertEqual(MyDetector.det3.zenith, 60.0)
    
    def test_from_full_name(self):
        self.assertEqual(MyDetector.from_full_name('Det0').name, 'det0')
        
        with self.assertRaises(ValueError):
            MyDetector.from_full_name('test')
    
    def test_from_num(self):
        self.assertEqual(MyDetector.from_num(0).name, 'det0')

        with self.assertRaises(ValueError):
            MyDetector.from_num(5)
    
    def test_from_str(self):
        self.assertEqual(MyDetector.from_str('det0').name, 'det0')

        with self.assertRaises(ValueError):
            MyDetector.from_str('test')

    def test_pointing(self):
        self.assertTupleEqual(MyDetector.det0.pointing(), (0.0, 15.0))
        
        
                
           
if __name__ == '__main__':
    unittest.main()


