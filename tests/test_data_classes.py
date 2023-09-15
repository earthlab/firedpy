import os
import shutil
import unittest
from datetime import datetime

from firedpy.data_classes import Base, BurnData

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))


class TestBase(unittest.TestCase):
    _test_data_dir = os.path.join(PROJECT_DIR, 'tests', 'test_data_dir')

    def tearDown(self) -> None:
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


if __name__ == '__main__':
    unittest.main()
