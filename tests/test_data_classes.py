import os
import re
import shutil
import unittest
from datetime import datetime
from unittest.mock import patch, Mock, MagicMock

import paramiko

from firedpy.data_classes import Base, BurnData

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))


class TestBase(unittest.TestCase):
    _test_data_dir = os.path.join(PROJECT_DIR, 'tests', 'test_data_dir')

    def tearDown(self) -> None:
        if os.path.exists(self._test_data_dir):
            shutil.rmtree(self._test_data_dir)

    def test_init(self):
        self.assertFalse(os.path.exists(self._test_data_dir))

        base = Base(out_dir=self._test_data_dir)

        self.assertIsNotNone(self._test_data_dir)
        self.assertTrue(os.path.exists(self._test_data_dir))
        self.assertIsInstance(base._date, str)
        self.assertIsInstance(datetime.strptime(base._date, '%m-%d-%Y'), datetime)
        self.assertTrue(base._cpus > 0)

        self.assertIsInstance(base._burn_area_dir, str)
        self.assertIsInstance(base._land_cover_dir, str)
        self.assertIsInstance(base._eco_region_raster_dir, str)
        self.assertIsInstance(base._eco_region_shapefile_dir, str)
        self.assertIsInstance(base._tables_dir, str)
        self.assertIsInstance(base._mosaics_dir, str)
        self.assertIsInstance(base._nc_dir, str)
        self.assertIsInstance(base._hdf_dir, str)

        self.assertIsInstance(base.MODIS_CRS, str)
        self.assertIsInstance(base.MODIS_SINUSOIDAL_PATH, str)
        self.assertIsInstance(base.CONUS_SHAPEFILE_PATH, str)
        self.assertIsInstance(base.DEFAULT_TILES, list)

    def test_initialize_save_dirs(self):
        self.assertFalse(any([
            os.path.exists(save_dir) for save_dir in [
                os.path.join(self._test_data_dir, "raster"),
                os.path.join(self._test_data_dir, "shape_files"),
                os.path.join(self._test_data_dir, 'rasters', 'burn_area'),
                os.path.join(self._test_data_dir, 'rasters', 'land_cover'),
                os.path.join(self._test_data_dir, 'rasters', 'eco_region'),
                os.path.join(self._test_data_dir, 'rasters', 'eco_region'),
                os.path.join(self._test_data_dir, 'tables'),
                os.path.join(self._test_data_dir, 'rasters', 'land_cover', 'mosaics'),
                os.path.join(self._test_data_dir, 'rasters', 'burn_area', 'netcdfs'),
                os.path.join(self._test_data_dir, 'rasters', 'burn_area', 'hdfs')
            ]
        ]))

        _ = Base(out_dir=self._test_data_dir)

        self.assertTrue(all([
            os.path.exists(save_dir) for save_dir in [
                os.path.join(self._test_data_dir, "rasters"),
                os.path.join(self._test_data_dir, "shape_files"),
                os.path.join(self._test_data_dir, 'rasters', 'burn_area'),
                os.path.join(self._test_data_dir, 'rasters', 'land_cover'),
                os.path.join(self._test_data_dir, 'rasters', 'eco_region'),
                os.path.join(self._test_data_dir, 'rasters', 'eco_region'),
                os.path.join(self._test_data_dir, 'tables'),
                os.path.join(self._test_data_dir, 'rasters', 'land_cover', 'mosaics'),
                os.path.join(self._test_data_dir, 'rasters', 'burn_area', 'netcdfs'),
                os.path.join(self._test_data_dir, 'rasters', 'burn_area', 'hdfs')
            ]
        ]))

    def test_get_shape_files(self):
        self.assertFalse(os.path.exists(os.path.join(self._test_data_dir, "shape_files",
                                                     'modis_sinusoidal_grid_world.shp')))
        self.assertFalse(os.path.exists(os.path.join(self._test_data_dir, "shape_files",
                                                     'conus.shp')))

        _ = Base(out_dir=self._test_data_dir)

        self.assertTrue(os.path.exists(os.path.join(self._test_data_dir, "shape_files",
                                                    'modis_sinusoidal_grid_world.shp')))
        self.assertTrue(os.path.exists(os.path.join(self._test_data_dir, "shape_files",
                                                    'conus.shp')))

    def test_convert_julian_date(self):
        # For January 1, 1970
        days = Base._convert_julian_date(1, 1970)
        self.assertEqual(days, 0)  # Expected 0 days since it's the base date

        # For January 2, 1970
        days = Base._convert_julian_date(2, 1970)
        self.assertEqual(days, 1)  # Expected 1 day since January 1, 1970

        # For January 1, 1971
        days = Base._convert_julian_date(1, 1971)
        self.assertEqual(days, 365)  # Expected 365 days (since 1970 was not a leap year)

        # For a random date like April 3, 1980
        date_1980 = datetime(1980, 4, 3)
        base_date = datetime(1970, 1, 1)
        expected_days = (date_1980 - base_date).days
        days = Base._convert_julian_date(94, 1980)
        self.assertEqual(days, expected_days)


class TestBurnData(unittest.TestCase):
    _test_data_dir = os.path.join(PROJECT_DIR, 'tests', 'test_data_dir')

    def tearDown(self) -> None:
        if os.path.exists(self._test_data_dir):
            shutil.rmtree(self._test_data_dir)

    def setUp(self) -> None:
        self.burn_data = BurnData(self._test_data_dir)
        self.hdfs = ["file1.hdf", "file2.hdf", "file3.hdf"]
        self.tile = "tile1"

    def test_init(self):
        self.assertIsNotNone(self.burn_data._base_sftp_folder)
        self.assertIsNotNone(self.burn_data._modis_template_path)
        self.assertTrue(self.burn_data._modis_template_path.startswith(self._test_data_dir))
        self.assertIsNotNone(self.burn_data._record_start_year)
        self.assertEqual(2000, self.burn_data._record_start_year)
        self.assertIsNotNone(self.burn_data._hdf_regex)

    def test_hdf_regex(self):
        burn = BurnData(out_dir=self._test_data_dir)

        nominal_file = 'MCD64A1.A2003032.h00v08.061.2021308094049.hdf'
        match = re.match(burn._hdf_regex, nominal_file)
        self.assertIsNotNone(match)
        group_dict = match.groupdict()
        self.assertTrue(
            all([param in group_dict.keys() for param in [
                'year', 'julian_day', 'horizontal_tile', 'vertical_tile', 'prod_year', 'prod_julian_day',
                'prod_hourminute', 'prod_second'
            ]])
        )
        self.assertEqual('2003', group_dict['year'])
        self.assertEqual('032', group_dict['julian_day'])
        self.assertEqual('00', group_dict['horizontal_tile'])
        self.assertEqual('08', group_dict['vertical_tile'])
        self.assertEqual('2021', group_dict['prod_year'])
        self.assertEqual('308', group_dict['prod_julian_day'])
        self.assertEqual('0940', group_dict['prod_hourminute'])
        self.assertEqual('49', group_dict['prod_second'])

        off_nominals = [
            'MCD64A1.A2003032.h00v08.061.2021308094049.nc4',
            'MCD64A1.A2003032.h005v08.061.2021308094049.hdf',
            'MCD12A5.A2003032.h00v08.061.2021308094049.hdf',
            'MCD64A1.A200032.h00v08.061.2021308094049.hdf',
            'MCD64A1.A2003032.h00v08.061.202130809404956.hdf',
            'MCD64A1.A2003032.h00v08.006.2021308094049.hdf'
        ]
        for off_nominal in off_nominals:
            self.assertIsNone(re.match(burn._hdf_regex, off_nominal))

    @patch('firedpy.data_classes.gdal.Open')  # Mock the gdal.Open method
    def test_verify_hdf_file_success(self, mock_open):
        mock_ds = Mock()  # Create a mock dataset
        mock_ds.GetSubDatasets.return_value = [1, 2, 3]  # Mock subdatasets
        mock_open.return_value = mock_ds  # Return the mock dataset when gdal.Open is called

        result = BurnData._verify_hdf_file('dummy_path.hdf')
        self.assertTrue(result)

    @patch('firedpy.data_classes.gdal.Open')
    def test_verify_hdf_file_failure_open(self, mock_open):
        mock_open.return_value = None

        result = BurnData._verify_hdf_file('dummy_path.hdf')
        self.assertFalse(result)

    @patch('firedpy.data_classes.gdal.Open')
    def test_verify_hdf_file_failure_subdatasets(self, mock_open):
        mock_ds = Mock()
        mock_ds.GetSubDatasets.return_value = []  # Empty subdatasets
        mock_open.return_value = mock_ds

        result = BurnData._verify_hdf_file('dummy_path.hdf')
        self.assertFalse(result)

    def test_generate_local_hdf_path(self):
        result = self.burn_data._generate_local_hdf_path("tile1", "hdf1.hdf")
        expected = os.path.join(self.burn_data._hdf_dir, 'tile1', 'hdf1.hdf')
        self.assertEqual(result, expected)

    def test_generate_remote_hdf_path(self):
        result = self.burn_data._generate_remote_hdf_path("tile1", "hdf1.hdf")
        expected = os.path.join('data', 'MODIS', 'C61', 'MCD64A1', 'HDF', 'tile1', 'hdf1.hdf')
        self.assertEqual(result, expected)

    def test_generate_remote_hdf_dir(self):
        result = self.burn_data._generate_remote_hdf_dir("tile1")
        expected = os.path.join('data', 'MODIS', 'C61', 'MCD64A1', 'HDF', 'tile1')
        self.assertEqual(result, expected)

    def test_generate_local_nc_path(self):
        result = self.burn_data._generate_local_nc_path("tile1")
        expected = os.path.join(self.burn_data._nc_dir, 'tile1.nc')
        self.assertEqual(result, expected)

    @patch("os.path.exists", return_value=True)
    @patch.object(BurnData, "_verify_hdf_file", return_value=True)
    def test_all_files_successful(self, mock_verify, mock_exists):
        mock_sftp_client = MagicMock()
        mock_sftp_client.get = MagicMock()

        self.burn_data._download_files(mock_sftp_client, self.tile, self.hdfs)
        self.assertEqual(mock_sftp_client.get.call_count, len(self.hdfs))

    @patch("os.path.exists", return_value=True)
    @patch.object(BurnData, "_verify_hdf_file", return_value=True)
    def test_some_files_fail_download(self, mock_verify, mock_exists):
        mock_sftp_client = MagicMock()
        mock_sftp_client.get = MagicMock()
        mock_sftp_client.get.side_effect = [Exception, None, None, None]

        self.burn_data._download_files(mock_sftp_client, self.tile, self.hdfs)
        self.assertEqual(mock_sftp_client.get.call_count, len(self.hdfs) + 1)  # One additional call due to one retry

    @patch.object(paramiko.SFTPClient, 'get')
    @patch("os.path.exists", return_value=True)
    @patch.object(BurnData, "_verify_hdf_file", side_effect=[False, True, True, True])
    def test_some_files_fail_verify(self, mock_verify, mock_exists, mock_sftp_get):
        mock_sftp_client = MagicMock()

        self.burn_data._download_files(mock_sftp_client, self.tile, self.hdfs)
        self.assertEqual(mock_verify.call_count, len(self.hdfs) + 1)  # One additional call due to one retry

    def test_never_succeed(self):
        mock_sftp_client = MagicMock()
        mock_sftp_client.get = MagicMock()
        mock_sftp_client.get.side_effect = Exception("Failed download")

        with self.assertRaises(IOError) as context:
            self.burn_data._download_files(mock_sftp_client, self.tile, self.hdfs)

        self.assertTrue(f'Error downloading burn data: max retries exceeded (3). Files not downloaded: {self.hdfs}'
                        in str(context.exception))


if __name__ == '__main__':
    unittest.main()
