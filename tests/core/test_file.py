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
import shutil
import unittest
import numpy as np
from tempfile import TemporaryDirectory, mkdtemp
from pathlib import Path
from astropy.io import fits
from gdt.core.file import FileContextManager, FitsFileContextManager


class TestFileContextManager(unittest.TestCase):

    def test_open(self):
        expected_contents = 'This is a test file.'

        with TemporaryDirectory() as temp_dir:
            test_file = Path(temp_dir) / 'test_file.txt'
            with test_file.open('w') as fp:
                fp.write(expected_contents)

            # Check to see if the file object is open within the context block
            with FileContextManager(test_file) as f:
                text = f.file_obj.read()
                self.assertEqual(expected_contents, text)
                self.assertFalse(f.file_obj.closed)

            # check to see if the file object is closed outside of context block
            self.assertTrue(f.file_obj.closed)


class TestFitsFileContextManager(unittest.TestCase):
    temp_dir: str
    fits_file: Path

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()

        table = fits.BinTableHDU.from_columns([
            fits.Column(name='SCLK_UTC', array=np.array([647049485.540078, 647049486.540079, 647049487.540079,
                                                         647049488.540078, 647049489.540078]), format='1D', unit='s'),
            fits.Column(name='QSJ_1', array=np.array([0.10640048601573592, 0.10659471756001844, 0.10678909044875709,
                                                      0.106984306723958, 0.10717798356382302]), format='1D'),
            fits.Column(name='QSJ_2', array=np.array([-0.06596779828740469, -0.06563509119982111, -0.0653022514085137,
                                                      -0.0649687052723567, -0.0646354981289785]), format='1D'),
            fits.Column(name='QSJ_3', array=np.array([-0.773998125107618, -0.7738618009494027, -0.7737249160835237,
                                                      -0.773587679290406, -0.7734507570667396]), format='1D'),
            fits.Column(name='QSJ_4', array=np.array([-0.620688398872268, -0.6208603015355695, -0.6210325759198554,
                                                      -0.621204900093449, -0.6213767445066635]), format='1D'),
            fits.Column(name='WSJ_1', array=np.array([-0.0007016791496425867, -0.0006979911704547703,
                                                      -0.0006990408292040229, -0.0006969416281208396,
                                                      -0.0006959060556255281]), format='1D', unit='rad/s'),
            fits.Column(name='WSJ_2', array=np.array([-0.00010338844731450081, -0.00010706065950216725,
                                                      -0.0001039113849401474, -0.00010601068788673729,
                                                      -0.00010811209358507767]), format='1D', unit='rad/s'),
            fits.Column(name='WSJ_3', array=np.array([-0.000534787483047694, -0.0005300534539856017,
                                                      -0.0005258548189885914, -0.0005300535704009235,
                                                      -0.0005353122833184898]), format='1D', unit='rad/s'),
            fits.Column(name='POS_X', array=np.array([-2165809.5, -2172932.0, -2180052.25, -2187169.75, -2194284.25]),
                        format='1E', unit='m'),
            fits.Column(name='POS_Y', array=np.array([5828483.0, 5825796.0, 5823102.5, 5820401.5, 5817693.5]),
                        format='1E', unit='m'),
            fits.Column(name='POS_Z', array=np.array([2979874.25, 2979937.25, 2979996.75, 2980052.5, 2980104.75]),
                        format='1E', unit='m'),
            fits.Column(name='VEL_X', array=np.array([-7123.95751953125, -7121.31982421875, -7118.67431640625,
                                                      -7116.0185546875, -7113.35400390625]), format='1E', unit='m/s'),
            fits.Column(name='VEL_Y', array=np.array([-2683.328857421875, -2690.4150390625, -2697.497314453125,
                                                      -2704.57666015625, -2711.653076171875]), format='1E', unit='m/s'),
            fits.Column(name='VEL_Z', array=np.array([64.80510711669922, 61.17137908935547, 57.5373649597168,
                                                      53.903621673583984, 50.27000045776367]), format='1E', unit='m/s'),
            fits.Column(name='SC_LAT', array=np.array([25.60576057434082, 25.606351852416992, 25.60690689086914,
                                                       25.607431411743164, 25.607921600341797]),
                        format='1E', unit='deg'),
            fits.Column(name='SC_LON', array=np.array([188.6570587158203, 188.7230224609375, 188.7889862060547,
                                                       188.85496520996094, 188.92092895507812]),
                        format='1E', unit='deg'),
            fits.Column(name='SADA_PY', array=np.array([-37.106075286865234, -37.11210250854492, -37.11532211303711,
                                                        -37.17781448364258, -37.29121017456055]),
                        format='1E', unit='deg'),
            fits.Column(name='SADA_NY', array=np.array([17.15755844116211, 17.227508544921875, 17.264917373657227,
                                                        17.2239933013916, 17.18284797668457]), format='1E', unit='deg'),
            fits.Column(name='FLAGS', array=np.array([1, 1, 1, 1, 1]), format='1I', bscale=1, bzero=32768)
        ])

        cls.temp_dir = mkdtemp()
        cls.fits_file = Path(cls.temp_dir) / 'test_data.fits'
        table.writeto(cls.fits_file)

    @classmethod
    def tearDownClass(cls) -> None:
        super().tearDownClass()
        shutil.rmtree(cls.temp_dir)

    def test_hdus(self):
        with FitsFileContextManager.open(self.fits_file) as f:
            self.assertEqual(f.num_hdus, 2)
            self.assertIsInstance(f.hdulist, fits.hdu.hdulist.HDUList)

    def test_filename(self):
        with FitsFileContextManager.open(self.fits_file) as f:
            self.assertEqual(f.filename, self.fits_file.name)

    def test_get_column_names(self):
        col_names = ('SCLK_UTC', 'QSJ_1', 'QSJ_2', 'QSJ_3', 'QSJ_4', 'WSJ_1',
                     'WSJ_2', 'WSJ_3', 'POS_X', 'POS_Y', 'POS_Z', 'VEL_X',
                     'VEL_Y', 'VEL_Z', 'SC_LAT', 'SC_LON', 'SADA_PY',
                     'SADA_NY', 'FLAGS')

        with FitsFileContextManager.open(self.fits_file) as f:
            self.assertTupleEqual(f.get_column_names(1), col_names)

    def test_column(self):
        lons = [188.6570587158203, 188.7230224609375, 188.7889862060547, 188.85496520996094, 188.92092895507812]
        with FitsFileContextManager.open(self.fits_file) as f:
            lon_list = f.column(1, 'SC_LON')[:5].tolist()
            for i in range(5):
                self.assertAlmostEqual(lon_list[i], lons[i], places=5)

    def test_columns_as_array(self):
        lons = [188.6570587158203, 188.7230224609375, 188.7889862060547, 188.85496520996094, 188.92092895507812]
        lats = [25.60576057434082, 25.606351852416992, 25.60690689086914, 25.607431411743164, 25.607921600341797]
        with FitsFileContextManager.open(self.fits_file) as f:
            scpos = f.columns_as_array(1, ['SC_LON', 'SC_LAT'])
            for i in range(5):
                self.assertAlmostEqual(scpos[i, 0], lons[i], places=5)
                self.assertAlmostEqual(scpos[i, 1], lats[i], places=5)

            scpos = f.columns_as_array(1, ['SC_LON', 'SC_LAT'], dtype=int)
            for i in range(5):
                self.assertEqual(scpos[i, 0], int(lons[i]))
                self.assertEqual(scpos[i, 1], int(lats[i]))

    def test_build_hdulist(self):
        with FitsFileContextManager.open(self.fits_file) as f:
            with self.assertRaises(NotImplementedError):
                f._build_hdulist()

    def test_write(self):
        with TemporaryDirectory() as this_path:
            with FitsFileContextManager.open(self.fits_file) as f:
                f.write(this_path)

            fpath = Path(this_path, self.fits_file.name)
            self.assertTrue(fpath.exists())

    def test_fsspec_cache_open_fits(self):
        trigger_name = "bn190915240"
        cache_path = self.temp_dir
        fsspec_kwargs = {"simplecache": {"cache_storage": cache_path, "same_names": True}}
        for fits_file in ["tcat", "trigdat"]:
            file_name = f"glg_{fits_file}_all_{trigger_name}_v0[0-9].fit"
            fits_url = f"simplecache::https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/20{trigger_name[2:4]}/{trigger_name}/current/{file_name}"
            with FitsFileContextManager.open(fits_url, use_fsspec=True, fsspec_kwargs=fsspec_kwargs) as f:
                assert f is not None
                assert f.num_hdus >= 1
                assert f.hdulist[0].name
                assert f.hdulist[0].header
                assert len(f.hdulist[0].header.cards) >= 1
                assert next(Path(cache_path).glob(file_name))
