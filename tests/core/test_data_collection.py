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
import numpy as np
import os
import unittest
from unittest import TestCase, skip
from gdt.core.collection import DataCollection

class DataObject():
    def __init__(self, filename, attribute):
        self.filename = filename
        self._attribute = attribute
    
    @property
    def attribute(self):
        return self._attribute
    
    def get_attribute(self):
        return self._attribute
    
    def set_attribute(self, attribute):
        self._attribute = attribute
    

class TestDataCollection(TestCase):
    d1 = DataObject('d1.dat', 1.0)
    d2 = DataObject('d2.dat', 2.0)
    d3 = DataObject('d3.dat', 3.0)
    collection = DataCollection.from_list([d1, d2, d3])
    
    def test_attributes(self):
        self.assertEqual(len(self.collection), 3)
        self.assertCountEqual(self.collection.items, [self.d1.filename, 
                                                      self.d2.filename,
                                                      self.d3.filename])
        self.assertEqual(self.collection.types, DataObject)

    def test_get_item(self):
        item = self.collection.get_item(self.collection.items[0])
        self.assertEqual(item, self.d1)
    
    def test_remove_and_include(self):
        # remove
        self.collection.remove(self.collection.items[2])
        self.assertEqual(len(self.collection), 2)
        self.assertCountEqual(self.collection.items, [self.d1.filename, 
                                                      self.d2.filename])
        # include
        self.collection.include(self.d3)
        self.test_attributes()
    
    def test_to_list(self):
        thelist = self.collection.to_list()
        self.assertCountEqual(thelist, [self.d1, self.d2, self.d3])
    
    def test_item_attributes(self):
        attributes = self.collection.attribute()
        self.assertCountEqual(attributes, [1.0, 2.0, 3.0])

    def test_item_methods(self):
        attributes = self.collection.get_attribute()
        self.assertCountEqual(attributes, [1.0, 2.0, 3.0])

        
if __name__ == '__main__':
    unittest.main()