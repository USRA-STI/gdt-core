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
import os
import unittest
from gdt.core.heasarc import FtpFinder, BrowseCatalog

this_dir = os.path.dirname(__file__)

# classes for unit tests

class MyFtpFinder(FtpFinder):
    _root = '/fermi/data/gbm/triggers'
    def _construct_path(self, str_trigger_num):
        year = '20' + str_trigger_num[0:2]
        path = os.path.join(self._root, year, 'bn' + str_trigger_num,
                            'current')
        return path
        
class MyCatalog(BrowseCatalog):
    def __init__(self, cache_path=this_dir, **kwargs):
        super().__init__(cache_path, table='batsegrb', **kwargs)

class MyBadCatalog(BrowseCatalog):
    def __init__(self, cache_path=this_dir, **kwargs):
        super().__init__(cache_path, table='badcat', **kwargs)

# ------------------------------------------------------------------------------

@unittest.skipIf(
    os.environ.get('SKIP_HEASARC_FTP_TESTS', False), 'Skipping HEASARC FTP tests'
)
class TestFtpFinder(unittest.TestCase):
    
    def tearDown(self):
        try:
            os.remove(os.path.join(this_dir, 
                                   'glg_trigdat_all_bn170817529_v01.fit'))
        except:
            pass
    
    def test_initialize(self):
        finder = MyFtpFinder('170817529')
        self.assertGreater(finder.num_files, 1)
        self.assertGreater(len(finder.files), 1)
        
        files = finder.filter('trigdat', 'fit')
        self.assertListEqual(files, ['glg_trigdat_all_bn170817529_v01.fit'])
        finder.get(this_dir, files)
        self.assertTrue(os.path.exists(os.path.join(this_dir, files[0])))
    
    def test_cd(self):
        finder = MyFtpFinder()
        finder.cd('170817529')
        self.assertGreater(finder.num_files, 1)
        self.assertGreater(len(finder.files), 1)

        files = finder.filter('trigdat', 'fit')
        finder.get(this_dir, files, verbose=False)
    
    def test_ls(self):
        finder = MyFtpFinder('170817529')
        self.assertListEqual(finder.files, finder.ls('170817529'))

    def test_errors(self):
        with self.assertRaises(ValueError):
            MyFtpFinder('oops i did it again')
            
        with self.assertRaises(FileNotFoundError):
            MyFtpFinder().ls('...and again')


class TestBrowseCatalog(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.catalog = MyCatalog(verbose=True)
    
    @classmethod
    def tearDownClass(cls):
        os.remove(os.path.join(this_dir, 'batsegrb.fit'))
        os.remove(os.path.join(this_dir, 'badcat.fit'))
    
    def test_num_rows(self):
        self.assertEqual(self.catalog.num_rows, 2702)
    
    def test_num_cols(self):
        self.assertEqual(self.catalog.num_cols, 47)

    def test_columns(self):
        self.assertTupleEqual(self.catalog.columns[:3], 
                             ('TRIGGER_NUM', 'NAME', 'RA'))
        self.assertTupleEqual(self.catalog.columns[-3:], 
                            ('THRESHOLD_1024', 'THRESHOLD_256', 'THRESHOLD_64'))
    
    def test_table(self):
        # get full table
        table = self.catalog.get_table()
        self.assertEqual(len(table.dtype), 47)
        self.assertEqual(table.size, 2702)

        # get one column
        table = self.catalog.get_table(columns=['RA'])
        self.assertTupleEqual(table.dtype.names, ('RA',))
        self.assertEqual(table.size, 2702)

        # get multiple column
        table = self.catalog.get_table(columns=['RA', 'TRIGGER_NUM'])
        self.assertTupleEqual(table.dtype.names, ('RA','TRIGGER_NUM'))
        self.assertEqual(table.size, 2702)
    
    def test_column_range(self):
        # range for an int
        self.assertTupleEqual(self.catalog.column_range('DAY_TRIGGER'), 
                              (8367, 11690))
        
        # range for a string
        # have to strip trailing spaces because HEASARC made a change that broke
        col_range = tuple([v.strip() for v in self.catalog.column_range('NAME')])
        self.assertTupleEqual(col_range, ('4B 910421', 'GRB 991229-'))

        # range for a float
        self.assertTupleEqual(self.catalog.column_range('SECONDS_TRIGGER'), 
                              (14.729, 86390.986))
    
    def test_slice(self):
    
        # slice for an int
        cat = self.catalog.slice('DAY_TRIGGER', lo=11690)
        self.assertEqual(cat.num_rows, 1)
        self.assertListEqual(cat.get_table(['DAY_TRIGGER'])['DAY_TRIGGER'].tolist(),
                             [11690])

        cat = self.catalog.slice('DAY_TRIGGER', hi=8367)
        self.assertEqual(cat.num_rows, 1)
        self.assertListEqual(cat.get_table(['DAY_TRIGGER'])['DAY_TRIGGER'].tolist(),
                             [8367])

        cat = self.catalog.slice('DAY_TRIGGER', lo=8367, hi=8369)
        self.assertEqual(cat.num_rows, 2)
        self.assertListEqual(cat.get_table(['DAY_TRIGGER'])['DAY_TRIGGER'].tolist(),
                             [8369, 8367])
        
        # slice for a string
        cat = self.catalog.slice('NAME', lo='GRB 991229')
        self.assertEqual(cat.num_rows, 2)
        # have to strip trailing spaces because HEASARC made a change that broke
        names = [name.strip() for name in cat.get_table(['NAME'])['NAME'].tolist()]
        self.assertListEqual(names, ['GRB 991229-', 'GRB 991229-'])
        
        cat = self.catalog.slice('NAME', hi='4B 910421')
        self.assertEqual(cat.num_rows, 1)
        # have to strip trailing spaces because HEASARC made a change that broke
        names = [name.strip() for name in cat.get_table(['NAME'])['NAME'].tolist()]
        self.assertListEqual(names, ['4B 910421'])

        cat = self.catalog.slice('NAME', lo='4B 910501', hi='4B 910505')
        self.assertEqual(cat.num_rows, 5)
        # have to strip trailing spaces because HEASARC made a change that broke
        names = [name.strip() for name in cat.get_table(['NAME'])['NAME'].tolist()]
        self.assertListEqual(names, ['4B 910502-', '4B 910502-', '4B 910501', 
                                     '4B 910503', '4B 910505'])

        # slice for a float
        cat = self.catalog.slice('RA', lo=359.0)
        for ra in cat.get_table(['RA'])['RA'].tolist():
            self.assertGreaterEqual(ra, 359.0)

        cat = self.catalog.slice('RA', hi=1.0)
        for ra in cat.get_table(['RA'])['RA'].tolist():
            self.assertLessEqual(ra, 1.0)

        cat = self.catalog.slice('RA', lo=5.0, hi=10.0)
        for ra in cat.get_table(['RA'])['RA'].tolist():
            self.assertTrue((ra >= 5.0) or (ra <= 10.0))
        
        # slice out of range - no data
        cat = self.catalog.slice('RA', hi=-1.0)
        self.assertEqual(cat.num_rows, 0)

    def test_slices(self):
        cat = self.catalog.slices([('RA', 5.0, 10.0), ('DEC', -90.0, -70.0)])
        table = cat.get_table(columns=['RA', 'DEC'])
        for i in range(table.size):
            self.assertTrue((table['RA'][i] >= 5.0) or (table['RA'][i] <= 10.0))
            self.assertTrue((table['DEC'][i] >= -90.0) or \
                            (table['DEC'][i] <= -70.0))

    def test_cached(self):
        cat = MyCatalog(cached=True)
        self.assertEqual(cat.num_rows, self.catalog.num_rows)        
        self.assertEqual(cat.num_cols, self.catalog.num_cols)
    
    def test_errors(self):
         with self.assertRaises(OSError):
            cat = MyBadCatalog()

if __name__ == '__main__':
    unittest.main()


