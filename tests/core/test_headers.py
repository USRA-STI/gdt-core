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
import warnings
from gdt.core.headers import Header, FileHeaders

# classes for unit tests

class MyHeader(Header):
    name = 'PRIMARY'
    keywords = [('STRING', 'hello', 'A defined string value'),
                ('INT', 1, 'An integer value'),
                ('FLOAT', 5.7, 'A float value'),
                ('BOOL', True, 'A Boolean value'),
                ('EXTNAME', '', 'Extension Name'),
                ('HY-PHEN', '', 'Hypehanated keyword')]

class MySecondHeader(Header):
    name = 'SECONDARY'
    keywords = [('ONE_KEY', '', 'A keyword'),
                ('DATE', '', 'The date'),
                ('COMMENT', 'blahblah', 'A comment')]

class MyEmptyHeader(Header):
    name = 'Empty'   
    
class MyBadHeader(Header):
    pass

class MyFileHeaders(FileHeaders):
    _header_templates = [MyHeader(), MySecondHeader()]

class MyBadFileHeaders(FileHeaders):
    pass

#-----------------------------------------------------------------------------

class TestHeader(unittest.TestCase):
    
    def setUp(self):
        self.header = MyHeader()
    
    def test_string(self):
        self.assertEqual(self.header['STRING'], 'hello')
        self.assertEqual(self.header.comments['STRING'], 
                         'A defined string value')
        
        # update value with another string
        self.header['STRING'] = 'Another string'
        self.assertEqual(self.header['STRING'], 'Another string')
        
        # update to value that is convertible to string
        self.header['STRING'] = 1000
        self.assertEqual(self.header['STRING'], '1000')
        
        # update to None
        self.header['STRING'] = None
        self.assertIsNone(self.header['STRING'])
    
    def test_int(self):
        self.assertEqual(self.header['INT'], 1)        
        self.assertEqual(self.header.comments['INT'], 
                         'An integer value')

        # update value with another int
        self.header['INT'] = 5
        self.assertEqual(self.header['INT'], 5)

        # update to value that is convertible to int
        self.header['INT'] = '1000'
        self.assertEqual(self.header['INT'], 1000)
        
        # update to value that is not convertible to int
        with self.assertRaises(TypeError):
            self.header['INT'] = 'hello'

    def test_float(self):
        self.assertEqual(self.header['FLOAT'], 5.7)        
        self.assertEqual(self.header.comments['FLOAT'], 
                         'A float value')

        # update value with another float
        self.header['FLOAT'] = 42.17
        self.assertEqual(self.header['FLOAT'], 42.17)

        # update to value that is convertible to float
        self.header['FLOAT'] = '1000'
        self.assertEqual(self.header['FLOAT'], 1000.0)
        
        # update to value that is not convertible to float
        with self.assertRaises(TypeError):
            self.header['FLOAT'] = 'hello'

    def test_bool(self):
        self.assertTrue(self.header['BOOL'])        
        self.assertEqual(self.header.comments['BOOL'], 
                         'A Boolean value')

        # update value with another boolean
        self.header['BOOL'] = False
        self.assertFalse(self.header['BOOL'])

        # update to value that is convertible to boolean
        self.header['BOOL'] = 1
        self.assertTrue(self.header['BOOL'])
    
    def test_access_error(self):
        with self.assertRaises(KeyError):
            self.header['MISSING']
    
    def test_add_new_keyword_error(self):
        with self.assertRaises(AttributeError):
            self.headers['NEW_KEY'] = 'test'

    def test_creator(self):
        creator = self.header.creator()
        self.assertEqual(creator[0], 'CREATOR')

    def test_empty(self):
        header = MyEmptyHeader()
        self.assertEqual(header.name, 'Empty')

    def test_assign_on_init(self):
        headers = MyHeader(string='goodbye', hy_phen='testing')
        self.assertEqual(headers['STRING'], 'goodbye')
        self.assertEqual(headers['HY-PHEN'], 'testing')

    def test_init_errors(self):
        with self.assertRaises(AttributeError):
            MyBadHeader()


class TestFileHeaders(unittest.TestCase):
    
    def setUp(self):
        self.headers = MyFileHeaders()
    
    def test_access_by_index(self):
        self.assertEqual(self.headers[0].name, 'PRIMARY')
        self.assertEqual(self.headers[1].name, 'SECONDARY')
        with self.assertRaises(IndexError):
            self.headers[5]

    def test_access_by_keyword(self):
        self.assertEqual(self.headers['PRIMARY'].name, 'PRIMARY')
        self.assertEqual(self.headers['SECONDARY'].name, 'SECONDARY')
        with self.assertRaises(KeyError):
            self.headers['TERTIARY']

    def test_num_headers(self):
        self.assertEqual(self.headers.num_headers, 2)
    
    def test_copy(self):
        headers = self.headers.copy()
        self.assertEqual(headers[0]['STRING'], self.headers[0]['STRING'])
        self.assertEqual(headers[1]['ONE_KEY'], self.headers[1]['ONE_KEY'])
    
    def test_keys(self):
        self.assertListEqual(self.headers.keys(), ['PRIMARY', 'SECONDARY'])

    def test_from_headers(self):
        headers = MyFileHeaders.from_headers([MyHeader(), MySecondHeader()])
        self.assertEqual(headers[0]['STRING'], self.headers[0]['STRING'])
        self.assertEqual(headers[1]['ONE_KEY'], self.headers[1]['ONE_KEY'])
        
        # wrong number of headers
        with self.assertRaises(ValueError):
             MyFileHeaders.from_headers([MyHeader()])
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            MyFileHeaders.from_headers([MyHeader(), MyHeader()])
            assert len(w) == 3

    def test_creator(self):
        creator = self.headers.creator()
        self.assertEqual(creator[0], 'CREATOR')

    def test_init_errors(self):
        with self.assertRaises(AttributeError):
            MyBadFileHeaders()
           
if __name__ == '__main__':
    unittest.main()


