import os
import re
import shutil
import unittest
from datetime import datetime
from unittest.mock import patch, Mock, MagicMock, call

import paramiko
import numpy as np
from netCDF4 import Dataset
import geopandas as gpd

from firedpy.data_classes import Base, BurnData, EcoRegion
from osgeo import gdal
import pandas as pd

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

    @patch("osgeo.ogr.Open")
    @patch("osgeo.gdal.GetDriverByName")
    @patch("osgeo.osr.SpatialReference")
    def test_rasterize_vector_data(self, mock_spatial_ref, mock_get_driver, mock_open):
        # Mock objects and methods
        mock_layer = Mock()
        mock_src_data = Mock(GetLayer=Mock(return_value=mock_layer))
        mock_open.return_value = mock_src_data

        mock_band = Mock()
        mock_trgt = Mock(GetRasterBand=Mock(return_value=mock_band), SetGeoTransform=Mock())
        mock_get_driver.return_value = Mock(Create=Mock(return_value=mock_trgt))

        mock_spatial_ref_instance = Mock(ExportToWkt=Mock(return_value="crs_wkt"))
        mock_spatial_ref.return_value = mock_spatial_ref_instance

        # Actual parameters you would use to call the method
        src = os.path.join('shapefiles', 'ecoregion', 'NA_CEC_Eco_Level3.shp')
        dst = os.path.join('rasters', 'ecoregion', 'NA_CEC_Eco_Level3.tif')
        attribute = "US_L3CODE"
        resolution = 100
        crs = "WGS 84"
        extent = [0, 0, 1000, 1000]  # example extent
        all_touch = False
        na = -9999

        # Call the method
        Base._rasterize_vector_data(src, dst, attribute, resolution, crs, extent, all_touch, na)

        mock_open.assert_called_with(src)
        mock_get_driver.assert_called_with("GTiff")
        mock_get_driver.return_value.Create.assert_called_with(dst, 10, 11, 1, gdal.GDT_Float32)
        mock_trgt.SetGeoTransform.assert_called_with((0, resolution, 0, 1000, 0, -resolution))
        mock_spatial_ref_instance.ImportFromWkt.assert_called_with(crs)
        mock_band.SetNoDataValue.assert_called_with(na)


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

        self.assertEqual(min(output_file.variables['time'][:]), 11262)
        self.assertEqual(max(output_file.variables['time'][:]), 19509)

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
        self.assertTrue(all([130 <= n <= 131 for n in non_nulls[0]]))
        self.assertEqual(min(values[non_nulls]), 15240)  # Sept 23rd 2011
        self.assertEqual(max(values[non_nulls]), 15257)

        self.assertEqual(['spatial_ref', 'proj4', 'geo_transform', 'grid_mapping_name', 'false_easting',
                          'false_northing', 'longitude_of_central_meridian', 'longitude_of_prime_meridian',
                          'semi_major_axis', 'inverse_flattening'],
                         list(vars(output_file.variables['crs']).keys()))

        self.assertEqual(output_file.variables['crs'].spatial_ref,
                         '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs')
        self.assertEqual(output_file.variables['crs'].proj4,
                         '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs')
        tolerance = 1e-2  # or whatever tolerance you deem appropriate
        self.assertTrue(np.allclose(output_file.variables['crs'].geo_transform,
                                    np.array([-1.89031588e+07, 4.63312717e+02,
                                              0.00000000e+00, -1.11195052e+06,
                                              0.00000000e+00, -4.63312717e+02]),
                                    atol=tolerance))
        self.assertEqual(output_file.variables['crs'].grid_mapping_name, 'sinusoidal')
        self.assertEqual(output_file.variables['crs'].false_easting, 0.0)
        self.assertEqual(output_file.variables['crs'].false_northing, 0.0)
        self.assertEqual(output_file.variables['crs'].longitude_of_central_meridian, 0.0)
        self.assertEqual(output_file.variables['crs'].longitude_of_prime_meridian, 0.0)
        self.assertEqual(output_file.variables['crs'].semi_major_axis, 6371007.181)
        self.assertEqual(output_file.variables['crs'].inverse_flattening, 0.0)


class TestEcoRegion(unittest.TestCase):

    def test_init_method(self):
        out_dir = "sample_out_dir"

        # Creating an instance of EcoRegion
        eco_region_instance = EcoRegion(out_dir)

        # Assertions to check if the attributes are set correctly
        self.assertEqual(eco_region_instance._eco_region_ftp_url,
                         'ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip')

        self.assertIsNotNone(eco_region_instance._eco_region_ftp_url)
        self.assertTrue(eco_region_instance._eco_region_ftp_url.startswith(out_dir))
        self.assertIsNotNone(eco_region_instance._project_eco_region_path)
        self.assertTrue(eco_region_instance._project_eco_region_path.startswith(out_dir))
        self.assertIsNotNone(eco_region_instance._eco_region_path)
        self.assertTrue(eco_region_instance._eco_region_path.startswith(out_dir))
        self.assertIsNotNone(eco_region_instance._eco_region_csv_path)
        self.assertTrue(eco_region_instance._eco_region_csv_path.startswith(out_dir))
        self.assertIsNotNone(eco_region_instance._eco_region_raster_path)
        self.assertTrue(eco_region_instance._eco_region_raster_path.startswith(out_dir))
        self.assertIsNotNone(eco_region_instance._ref_cols)
        self.assertTrue(eco_region_instance._ref_cols.startswith(out_dir))

        # Check if eco_region_data_frame is initialized to None
        self.assertIsNone(eco_region_instance.eco_region_data_frame)

    def test_basic_capitalization(self):
        self.assertEqual(
            EcoRegion._normalize_string("hello world"),
            "Hello World"
        )

    def test_special_characters(self):
        self.assertEqual(
            EcoRegion._normalize_string("hello/world good-morning"),
            "Hello/World Good-Morning"
        )

    def test_special_word_and(self):
        self.assertEqual(
            EcoRegion._normalize_string("hello and world"),
            "Hello and World"
        )

    def test_special_word_usa(self):
        self.assertEqual(
            EcoRegion._normalize_string("hello USA world"),
            "Hello USA World"
        )

    def test_combination(self):
        self.assertEqual(
            EcoRegion._normalize_string("hello and USA/world good-morning"),
            "Hello and USA/World Good-Morning"
        )

    @patch("geopandas.read_file")
    @patch("os.path.exists")
    def test_read_from_ftp(self, mock_exists, mock_read_file):
        mock_exists.return_value = False
        mock_read_file.return_value = gpd.GeoDataFrame()

        eco_region = EcoRegion('test_out_dir')
        df = eco_region._read_eco_region_file()

        mock_read_file.assert_called_with(eco_region._eco_region_ftp_url)
        self.assertIsInstance(df, gpd.GeoDataFrame)

    @patch("geopandas.read_file")
    @patch("os.path.exists")
    @patch("shutil.copy")
    def test_read_from_local_when_ftp_fails(self, mock_copy, mock_exists, mock_read_file):
        mock_exists.return_value = False
        mock_read_file.side_effect = [Exception("FTP Error"), gpd.GeoDataFrame()]

        eco_region = EcoRegion('test_out_dir')
        df = eco_region._read_eco_region_file()

        mock_read_file.assert_called_with(eco_region._eco_region_path)
        self.assertIsInstance(df, gpd.GeoDataFrame)

    @patch("geopandas.read_file")
    @patch("os.path.exists")
    @patch("shutil.copy")
    def test_handle_total_failure(self, mock_copy, mock_exists, mock_read_file):
        mock_exists.return_value = False
        mock_read_file.side_effect = Exception("Total Failure")

        eco_region = EcoRegion('test_out_dir')
        with self.assertRaises(Exception) as context:
            eco_region._read_eco_region_file()

        self.assertEqual(str(context.exception), "Total Failure")

    @patch.object(EcoRegion, '_read_eco_region_file')
    @patch.object(EcoRegion, '_normalize_string')
    @patch("pandas.DataFrame.to_csv")
    def test_get_eco_region(self, mock_to_csv, mock_normalize, mock_read_file):
        # Mock the data returned by `_read_eco_region_file`
        mock_df = pd.DataFrame({
            'NA_L3CODE': ['A', 'B', 'C'],
            'NA_L3NAME': ['region_A', 'region_B', 'region_C'],
            # Add other columns as per `_ref_cols`
        })
        mock_read_file.return_value = mock_df

        # Mock `_normalize_string` method to just return the input string for simplicity
        mock_normalize.side_effect = lambda x: x

        # Initialize EcoRegion and call `get_eco_region`
        eco_region = EcoRegion('your_out_dir')
        eco_region.get_eco_region()

        # Check that `_read_eco_region_file` was called once
        mock_read_file.assert_called_once()

        # Check that `_normalize_string` was called for each element in `eco_ref`
        calls = [call(cell) for cell in mock_df[eco_region._ref_cols].values.flatten()]
        mock_normalize.assert_has_calls(calls, any_order=True)

        # Check that `to_csv` was called once with the correct path
        mock_to_csv.assert_called_once_with(eco_region._eco_region_csv_path, index=False)

        # Check that `eco_region_data_frame` is assigned correctly
        pd.testing.assert_frame_equal(eco_region.eco_region_data_frame, mock_df.applymap(str))

    @patch("geopandas.read_file")
    @patch.object(EcoRegion, '_rasterize_vector_data')
    @patch("os.path.exists", return_value=True)  # Mocking all `os.path.exists` calls to return True
    @patch("gdal.Open")
    @patch("os.path.join", return_value="mock_path")  # Mocking `os.path.join` to always return "mock_path"
    @patch("glob.glob",
           return_value=["mock_hdf_file_path"])  # Mocking `glob.glob` to always return a list with one mock path
    def test_create_eco_region_raster(self, mock_glob, mock_join, mock_gdal_open,
                                      mock_exists, mock_rasterize, mock_read_file):
        # Mocking data returned by `read_file`
        mock_df = gpd.GeoDataFrame({
            'h': [1, 2, 3],
            'v': [4, 5, 6],
            'tile': ['h01v04', 'h02v05', 'h03v06'],
        })
        mock_read_file.return_value = mock_df

        # Mocking data returned by `gdal.Open`
        mock_file_pointer = Mock()
        mock_ds = Mock()
        mock_file_pointer.GetSubDatasets.return_value = [(None, None)]
        mock_ds.GetGeoTransform.return_value = (0, 1, 0, 0, 0, -1)
        mock_ds.RasterXSize = 1
        mock_ds.RasterYSize = 1
        mock_ds.GetProjection.return_value = "mock_projection"
        mock_gdal_open.return_value = mock_ds

        # Setting up a mock for the EcoRegion instance's `_generate_local_burn_hdf_dir` method
        mock_eco_region = EcoRegion('test_out_dir')
        mock_eco_region._generate_local_burn_hdf_dir = Mock(return_value="mock_burn_dir")

        # Forcing `get_eco_region` to assign a non-None value to `eco_region_data_frame` to skip its call
        mock_eco_region.eco_region_data_frame = pd.DataFrame()

        # Calling the method with the mock objects and check behaviors
        tiles = ['h01v04', 'h03v06']
        mock_eco_region.create_eco_region_raster(tiles)

        # Check the `read_file` is called correctly
        mock_read_file.assert_called_once_with(mock_eco_region._modis_sinusoidal_grid_shape_path)

        # Check the `_generate_local_burn_hdf_dir` is called for each `extent_tile`
        mock_eco_region._generate_local_burn_hdf_dir.assert_has_calls(
            [call(tile) for tile in ['h01v04', 'h01v04', 'h03v06', 'h03v06']])

        # Check the `_rasterize_vector_data` is called correctly
        mock_rasterize.assert_called_once_with(mock_eco_region.eco_region_data_frame,
                                               mock_eco_region._eco_region_raster_path,
                                               "US_L3CODE", 1, "mock_projection", [0, -1, 1, 0])

    @patch("geopandas.read_file")
    @patch.object(EcoRegion, '_rasterize_vector_data')
    @patch.object(EcoRegion, 'get_eco_region')
    @patch("os.path.exists", return_value=False)  # Mocking all `os.path.exists` calls to return True
    @patch("gdal.Open")
    @patch("os.path.join", return_value="mock_path")  # Mocking `os.path.join` to always return "mock_path"
    @patch("glob.glob",
           return_value=[
               "mock_hdf_file_path"])  # Mocking `glob.glob` to always return a list with one mock path
    def test_create_eco_region_raster(self, mock_glob, mock_join, mock_gdal_open,
                                      mock_exists, mock_get_eco_region, mock_rasterize, mock_read_file):
        mock_df = gpd.GeoDataFrame({
            'h': [1, 2, 3],
            'v': [4, 5, 6],
            'tile': ['h01v04', 'h02v05', 'h03v06'],
        })

        mock_get_eco_region.return_value = mock_df

        # Mocking data returned by `gdal.Open`
        mock_file_pointer = Mock()
        mock_ds = Mock()
        mock_file_pointer.GetSubDatasets.return_value = [(None, None)]
        mock_ds.GetGeoTransform.return_value = (0, 1, 0, 0, 0, -1)
        mock_ds.RasterXSize = 1
        mock_ds.RasterYSize = 1
        mock_ds.GetProjection.return_value = "mock_projection"
        mock_gdal_open.return_value = mock_ds

        # Setting up a mock for the EcoRegion instance's `_generate_local_burn_hdf_dir` method
        mock_eco_region = EcoRegion('test_out_dir')
        mock_eco_region._generate_local_burn_hdf_dir = Mock(return_value="mock_burn_dir")
        mock_eco_region.eco_region_data_frame = None

        # Calling the method with the mock objects and check behaviors
        tiles = ['h01v04', 'h03v06']
        mock_eco_region.create_eco_region_raster(tiles)

        # Check the `read_file` is called correctly
        mock_read_file.assert_called_once_with(mock_eco_region._modis_sinusoidal_grid_shape_path)

        # Check the `_generate_local_burn_hdf_dir` is called for each `extent_tile`
        mock_eco_region._generate_local_burn_hdf_dir.assert_has_calls(
            [call(tile) for tile in ['h01v04', 'h01v04', 'h03v06', 'h03v06']])

        # Check the `_rasterize_vector_data` is called correctly
        mock_rasterize.assert_called_once_with(mock_eco_region.eco_region_data_frame,
                                               mock_eco_region._eco_region_raster_path,
                                               "US_L3CODE", 1, "mock_projection", [0, -1, 1, 0])

        self.assertTrue(mock_get_eco_region.call_count, 1)


if __name__ == '__main__':
    unittest.main()
