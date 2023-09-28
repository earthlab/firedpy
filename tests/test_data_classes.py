import os
import re
import shutil
import unittest
from datetime import datetime
from unittest.mock import patch, Mock, MagicMock

import paramiko
import numpy as np
from netCDF4 import Dataset

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

    def test_convert_ordinal_date(self):
        # For January 1, 1970
        days = Base._convert_ordinal_to_unix_day(1970, 1)
        self.assertEqual(days, 0)  # Expected 0 days since it's the base date

        # For January 2, 1970
        days = Base._convert_ordinal_to_unix_day(1970, 2)
        self.assertEqual(days, 1)  # Expected 1 day since January 1, 1970

        # For January 1, 1971
        days = Base._convert_ordinal_to_unix_day(1971, 1)
        self.assertEqual(days, 365)  # Expected 365 days (since 1970 was not a leap year)

        # For a random date like April 3, 1980
        date_1980 = datetime(1980, 4, 3)
        base_date = datetime(1970, 1, 1)
        expected_days = (date_1980 - base_date).days
        days = Base._convert_ordinal_to_unix_day(1980, 94)
        self.assertEqual(days, expected_days)

    def test_convert_dates(self):
        base_class = Base(self._test_data_dir)
        arr = np.array([[0, 0], [0, 365], [1, 0]])
        year = 1971
        expected = np.array([[0, 0], [0, 365 + 364], [365, 0]])
        result = base_class._convert_dates(arr, year)
        np.testing.assert_array_equal(result, expected)


class TestBurnData(unittest.TestCase):
    _test_data_dir = os.path.join(PROJECT_DIR, 'tests', 'test_data_dir')

    def tearDown(self) -> None:
        if os.path.exists(self._test_data_dir):
            shutil.rmtree(self._test_data_dir)

    def setUp(self) -> None:
        if os.path.exists(self._test_data_dir):
            shutil.rmtree(self._test_data_dir)
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

    @patch("os.path.exists", return_value=False)
    @patch.object(BurnData, "_verify_hdf_file", return_value=True)
    def test_all_files_successful(self, mock_verify, mock_exists):
        mock_sftp_client = MagicMock()
        mock_sftp_client.get = MagicMock()

        self.burn_data._download_files(mock_sftp_client, self.tile, self.hdfs)
        self.assertEqual(mock_sftp_client.get.call_count, len(self.hdfs))

    @patch("os.path.exists", return_value=False)
    @patch.object(BurnData, "_verify_hdf_file", return_value=True)
    def test_some_files_fail_download(self, mock_verify, mock_exists):
        mock_sftp_client = MagicMock()
        mock_sftp_client.get = MagicMock()
        mock_sftp_client.get.side_effect = [Exception, None, None, None]

        self.burn_data._download_files(mock_sftp_client, self.tile, self.hdfs)
        self.assertEqual(mock_sftp_client.get.call_count, len(self.hdfs) + 1)  # One additional call due to one retry

    @patch.object(paramiko.SFTPClient, 'get')
    @patch("os.path.exists", return_value=False)
    @patch.object(BurnData, "_verify_hdf_file", side_effect=[False, True, True, True])
    def test_some_files_fail_verify(self, mock_verify, mock_exists, mock_sftp_get):
        mock_sftp_client = MagicMock()

        self.burn_data._download_files(mock_sftp_client, self.tile, self.hdfs)
        self.assertEqual(mock_verify.call_count, len(self.hdfs) + 1)  # One additional call due to one retry

    @patch("os.path.exists", return_value=False)
    def test_never_succeed(self, mock_os):
        mock_sftp_client = MagicMock()
        mock_sftp_client.get = MagicMock()
        mock_sftp_client.get.side_effect = Exception("Failed download")

        with self.assertRaises(IOError) as context:
            self.burn_data._download_files(mock_sftp_client, self.tile, self.hdfs)

        print(str(context.exception))

        self.assertTrue(f'Error downloading burn data: max retries exceeded (3). Files not downloaded or not able to '
                        f'open: {self.hdfs}' in str(context.exception))

    def test_extract_date_parts_valid(self):
        filename = "MCD64A1.A2022123.h12v34.061.2023123123456.hdf"
        result = self.burn_data._extract_date_parts(filename)
        self.assertEqual(result, (2022, 123))

    def test_extract_date_parts_invalid(self):
        filename = "MCD64X1.A2022123.h12v34.061.20231231234567.hdf"  # Intentionally changed A to X
        result = self.burn_data._extract_date_parts(filename)
        self.assertIsNone(result)

    def test_write_nc_nominal(self):
        burn_data = BurnData(self._test_data_dir)
        tile = 'h01v10'
        shutil.copytree(os.path.join(PROJECT_DIR, 'tests', 'data', 'burn_data', tile),
                        burn_data._generate_local_hdf_dir(tile))

        output_nc = os.path.join(self._test_data_dir, burn_data._generate_local_nc_path(tile))
        self.assertFalse(os.path.exists(output_nc))
        burn_data._write_ncs([tile])
        self.assertTrue(os.path.exists(output_nc))

        output_file = Dataset(output_nc)
        self.assertTrue(
            all(s in vars(output_file).keys() for s in ['title', 'subtitle', 'description', 'date', 'projection',
                                                        'Conventions']))
        self.assertEqual(output_file.title, 'Burn Days')
        self.assertEqual(output_file.subtitle, 'Burn Days Detection by MODIS since 1970.')
        self.assertEqual(output_file.description, 'The day that a fire is detected.')
        self.assertEqual(output_file.date, datetime.today().strftime('%Y-%m-%d')
                         )
        self.assertEqual(output_file.projection, 'MODIS Sinusoidal')
        self.assertEqual(output_file.Conventions, 'CF-1.6')

        self.assertEqual(['y', 'x', 'time', 'value', 'crs'], list(output_file.variables.keys()))

        self.assertEqual(['standard_name', 'long_name', 'units'], list(vars(output_file.variables['y']).keys()))

        self.assertEqual(output_file.variables['y'].standard_name, 'projection_y_coordinate')
        self.assertEqual(output_file.variables['y'].long_name, 'y coordinate of projection')
        self.assertEqual(output_file.variables['y'].units, 'm')

        self.assertAlmostEquals(min(output_file.variables['y'][:]), -2223437.7266224716)
        self.assertAlmostEquals(max(output_file.variables['y'][:]), -1111950.519672)

        self.assertEqual(['standard_name', 'long_name', 'units'], list(vars(output_file.variables['x']).keys()))

        self.assertEqual(output_file.variables['x'].standard_name, 'projection_x_coordinate')
        self.assertEqual(output_file.variables['x'].long_name, 'x coordinate of projection')
        self.assertEqual(output_file.variables['x'].units, 'm')

        self.assertAlmostEquals(min(output_file.variables['x'][:]), -18903158.834333)
        self.assertAlmostEquals(max(output_file.variables['x'][:]), -17791671.627382528)

        self.assertEqual(['units', 'standard_name', 'calendar'], list(vars(output_file.variables['time']).keys()))

        self.assertEqual(output_file.variables['time'].standard_name, 'time')
        self.assertEqual(output_file.variables['time'].calendar, 'gregorian')
        self.assertEqual(output_file.variables['time'].units, 'days since 1970-01-01')

        self.assertEqual(min(output_file.variables['time'][:]), 18262)
        self.assertEqual(max(output_file.variables['time'][:]), 18597)

        self.assertEqual(['_FillValue', 'standard_name', 'long_name', 'grid_mapping'],
                         list(vars(output_file.variables['value']).keys()))

        self.assertEqual(output_file.variables['value'].standard_name, 'day')
        self.assertEqual(output_file.variables['value']._FillValue, -9999)
        self.assertEqual(output_file.variables['value'].long_name, 'Burn Days')
        self.assertEqual(output_file.variables['value'].grid_mapping, 'crs')

        self.assertEqual(output_file.variables['value'][:].shape, (272, 2400, 2400))
        values = output_file.variables['value'][:]
        non_nulls = np.where(values > 0)

        # All months between 130 and 131 for burn events
        self.assertTrue(all([130 >= n >= 131 for n in non_nulls[0]]))
        self.assertEqual(min(values[non_nulls]), 15240)  # Sept 23rd 2011
        self.assertEqual(min(values[non_nulls]), 15257)

        self.assertEqual(['spatial_ref', 'proj4', 'geo_transform', 'grid_mapping_name', 'false_easting',
                          'false_northing', 'longitude_of_central_meridian', 'longitude_of_prime_meridian',
                          'semi_major_axis', 'inverse_flattening'],
                         list(vars(output_file.variables['crs']).keys()))

        self.assertEqual(output_file.variables['crs'].spatial_ref,
                         '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs')
        self.assertEqual(output_file.variables['crs'].proj4,
                         '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs')
        self.assertEqual(output_file.variables['crs'].geo_transform, np.array([-1.89031588e+07, 4.63312717e+02,
                                                                               0.00000000e+00, -1.11195052e+06,
                                                                               0.00000000e+00, -4.63312717e+02]))
        self.assertEqual(output_file.variables['crs'].grid_mapping_name, 'sinusoidal')
        self.assertEqual(output_file.variables['crs'].false_easting, 0.0)
        self.assertEqual(output_file.variables['crs'].false_northing, 0.0)
        self.assertEqual(output_file.variables['crs'].longitude_of_central_meridian, 0.0)
        self.assertEqual(output_file.variables['crs'].longitude_of_prime_meridian, 0.0)
        self.assertEqual(output_file.variables['crs'].semi_major_axis, 6371007.181)
        self.assertEqual(output_file.variables['crs'].inverse_flattening, 0.0)


if __name__ == '__main__':
    unittest.main()
