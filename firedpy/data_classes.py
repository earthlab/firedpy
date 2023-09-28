import datetime as dt
import os
import re
import shutil
import sys
from datetime import datetime
from glob import glob
from typing import List, Union
import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import paramiko
import rasterio
from netCDF4 import Dataset
from osgeo import gdal, osr
from rasterio.merge import merge
from tqdm import tqdm

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))

logging.basicConfig(filename='app.log', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')


class Base:
    """
    """
    MODIS_CRS = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    MODIS_SINUSOIDAL_PATH = os.path.join(PROJECT_DIR, "ref", "modis_grid.gpkg")
    CONUS_SHAPEFILE_PATH = os.path.join(PROJECT_DIR, 'ref', 'boundaries', 'conus.gpkg')
    DEFAULT_TILES = ["h08v04", "h09v04", "h10v04", "h11v04", "h12v04",
                     "h13v04", "h08v05", "h09v05", "h10v05", "h11v05",
                     "h12v05", "h08v06", "h09v06", "h10v06", "h11v06"]

    def __init__(self, out_dir: str):
        os.makedirs(out_dir, exist_ok=True)

        self._out_dir = out_dir
        self._date = datetime.today().strftime('%m-%d-%Y')
        self._cpus = os.cpu_count()

        self._raster_dir = os.path.join(out_dir, 'rasters')
        self._shape_file_dir = os.path.join(self._out_dir, "shape_files")
        self._burn_area_dir = os.path.join(self._raster_dir, 'burn_area')
        self._land_cover_dir = os.path.join(self._raster_dir, 'land_cover')
        self._eco_region_raster_dir = os.path.join(self._raster_dir, 'eco_region')
        self._eco_region_shapefile_dir = os.path.join(self._raster_dir, 'eco_region')
        self._tables_dir = os.path.join(out_dir, 'tables')

        self._mosaics_dir = os.path.join(self._land_cover_dir, 'mosaics')
        self._nc_dir = os.path.join(self._burn_area_dir, 'netcdfs')
        self._hdf_dir = os.path.join(self._burn_area_dir, 'hdfs')

        # Initialize output directory folders and files
        self._initialize_save_dirs()
        self._get_shape_files()

    def _initialize_save_dirs(self):
        for f in [
            self._raster_dir,
            self._shape_file_dir,
            self._burn_area_dir,
            self._land_cover_dir,
            self._eco_region_raster_dir,
            self._eco_region_shapefile_dir,
            self._tables_dir,
            self._mosaics_dir,
            self._nc_dir,
            self._hdf_dir
        ]:
            os.makedirs(f, exist_ok=True)

    def _get_shape_files(self):
        """
        Just to grab some basic shapefiles needed for calculating statistics.
        """

        files_to_copy = {
            'modis_sinusoidal_grid_world.shp': self.MODIS_SINUSOIDAL_PATH,
            'conus.shp': self.CONUS_SHAPEFILE_PATH
        }

        for dest_name, source_path in files_to_copy.items():
            dest_path = os.path.join(self._shape_file_dir, dest_name)
            if not os.path.exists(dest_path):
                shutil.copy(source_path, dest_path)

    @staticmethod
    def _convert_julian_date(year: int, julian_day: int) -> int:
        base = dt.datetime(1970, 1, 1)
        date = dt.datetime(year, 1, 1) + dt.timedelta(int(julian_day - 1))
        days = date - base
        return days.days

    def _convert_dates(self, array, year) -> np.array:
        """Convert every day in an array to days since Jan 1 1970"""
        # Loop through each position with data and convert
        ys, xs = np.where(array > 0)

        for y, x in zip(ys, xs):
            array[y, x] = self._convert_julian_date(year, array[y, x])

        return array


class BurnData(Base):
    def __init__(self, out_dir: str):
        super().__init__(out_dir)
        self._base_sftp_folder = os.path.join('data', 'MODIS', 'C61', 'MCD64A1', 'HDF')
        self._modis_template_path = os.path.join(out_dir, 'rasters', 'mosaic_template.tif')
        self._record_start_year = 2000
        self._hdf_regex = r'MCD64A1\.A(?P<year>\d{4})(?P<julian_day>\d{3})\.h(?P<horizontal_tile>\d{2})v(?P<vertical_tile>\d{2})\.061\.(?P<prod_year>\d{4})(?P<prod_julian_day>\d{3})(?P<prod_hourminute>\d{4})(?P<prod_second>\d{2})\.hdf$'

    @staticmethod
    def _verify_hdf_file(file_path) -> bool:
        # Open the HDF file
        hdf_ds = gdal.Open(file_path)

        if hdf_ds is None:
            print(f"Failed to open {file_path}.")
            return False

        # List available sub-datasets (specific to HDF)
        sub_datasets = hdf_ds.GetSubDatasets()

        if not sub_datasets:
            print("No sub-datasets found.")
            return False

        return True

    def _generate_local_hdf_path(self, tile: str, hdf_name: str) -> str:
        return os.path.join(self._hdf_dir, tile, hdf_name)

    def _generate_local_hdf_dir(self, tile: str) -> str:
        return os.path.join(self._hdf_dir, tile)

    def _generate_remote_hdf_path(self, tile: str, hdf_name: str) -> str:
        return os.path.join(self._base_sftp_folder, tile, hdf_name)

    def _generate_remote_hdf_dir(self, tile: str) -> str:
        return os.path.join(self._base_sftp_folder, tile)

    def _generate_local_nc_path(self, tile: str) -> str:
        return os.path.join(self._nc_dir, f"{tile}.nc")

    def _download_files(self, sftp_client: paramiko.SFTPClient, tile: str, hdfs: List[str],
                        max_retries: int = 3) -> None:
        attempt = 0
        retries = hdfs.copy()
        while attempt < max_retries and retries:
            next_retries = []
            for hdf_file in tqdm(retries):
                remote_path = self._generate_remote_hdf_path(tile, hdf_file)
                local_path = self._generate_local_hdf_path(tile, hdf_file)
                os.makedirs(os.path.dirname(local_path), exist_ok=True)
                try:
                    sftp_client.get(remote_path, local_path)
                    if not self._verify_hdf_file(local_path):
                        next_retries.append(hdf_file)
                except Exception as e:
                    next_retries.append(hdf_file)

            retries = next_retries
            attempt += 1

        if attempt == max_retries and retries:
            print('Raising exception')
            raise IOError(f'Error downloading burn data: max retries exceeded ({max_retries}). Files not downloaded: {retries}')

    def get_burns(self, tiles: List[str], start_year: int = None, end_year: int = None):
        """
        This will download the MODIS burn event data set tiles and create a
        singular mosaic to use as a template file for coordinate reference
        information and geometries.

        User manual:
            http://modis-fire.umd.edu/files/MODIS_C6_BA_User_Guide_1.2.pdf

        Update 02/2021 -> fuoco server transitioned to SFTP Dec 2020
            Update firedpy to use Paramiko SSHClient / SFTPClient
            Server-side changes are described in the user manual linked above

        SFTP:
            sftp://fire:burnt@fuoco.geog.umd.edu/gfed4/MCD64A1/C6/
            username: fire
            password: burnt

        """
        # Check into the UMD SFTP fuoco server using Paramiko
        ssh_client = paramiko.SSHClient()
        ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh_client.connect(hostname="fuoco.geog.umd.edu", username="fire", password="burnt")
        print("Connected to 'fuoco.geog.umd.edu' ...")

        year_range = np.arange(start_year if start_year is not None else self._record_start_year,
                               end_year if end_year is not None else datetime.now().year + 1, 1)

        # Open the connection to the SFTP
        sftp_client = ssh_client.open_sftp()

        for tile in tiles:
            nc_file = self._generate_local_nc_path(tile)

            if os.path.exists(nc_file):
                continue

            print(f"Downloading/Checking HDF files for: {tile}")
            hdf_files = []
            for hdf_file in sftp_client.listdir(self._generate_remote_hdf_dir(tile)):
                match = re.match(self._hdf_regex, hdf_file)
                if match is not None and int(match.groupdict()['year']) in year_range:
                    hdf_files.append(hdf_file)

            if not hdf_files:
                print(f"No MCD64A1 Product for tile: {tile}, skipping...")

            try:
                self._download_files(sftp_client, tile, hdf_files)
            except IOError as e:
                continue

        # Close new SFTP connection
        ssh_client.close()
        sftp_client.close()
        print("Disconnected from 'fuoco.geog.umd.edu' ...")

        if not os.path.exists(self._modis_template_path):
            self._write_modis_template_file()

        self._write_ncs(tiles)

    def _write_modis_template_file(self):
        # Merge one year into a reference mosaic
        print("Creating reference mosaic ...")
        folders = glob(os.path.join(self._hdf_dir, "*"))
        file_groups = [glob(os.path.join(f, "*hdf")) for f in folders]
        for f in file_groups:
            f.sort()
        files = [f[0] for f in file_groups]
        dss = [rasterio.open(f).subdatasets[0] for f in files]
        tiles = [rasterio.open(d) for d in dss]
        mosaic, transform = merge(tiles)
        crs = tiles[0].meta.copy()
        crs.update({"driver": "GTIFF",
                    "height": mosaic.shape[1],
                    "width": mosaic.shape[2],
                    "transform": transform})

        with rasterio.open(self._modis_template_path, "w+", **crs) as dst:
            dst.write(mosaic)

    def _extract_date_parts(self, filename):
        filename = os.path.basename(filename)
        match = re.match(self._hdf_regex, filename)
        if match:
            return int(match.groupdict()['year']), int(match.groupdict()['julian_day'])
        return None

    def _write_ncs(self, tiles: List[str]):
        """
        Take in a time series of files for the MODIS burn detection dataset and
        create a singular netcdf file.
        """
        # Build the net cdfs here
        for tile_id in tiles:
            try:
                hdf_dir = self._generate_local_hdf_dir(tile_id)

                files = sorted([os.path.join(hdf_dir, f) for f in os.listdir(hdf_dir) if self._extract_date_parts(f)
                                is not None], key=self._extract_date_parts)

                if not files:
                    print(f'No hdf files for tile {tile_id} in {hdf_dir}')

                nc_file_name = self._generate_local_nc_path(tile_id)

                # Skip if it exists already
                if os.path.exists(nc_file_name):
                    print(tile_id + " netCDF file exists, skipping...")
                else:
                    # Use a sample to get geography information and geometries
                    print("Building netcdf for tile " + tile_id)
                    sample = files[0]
                    ds = gdal.Open(sample).GetSubDatasets()[0][0]
                    hdf = gdal.Open(ds)
                    geom = hdf.GetGeoTransform()
                    proj = hdf.GetProjection()
                    data = hdf.GetRasterBand(1)
                    crs = osr.SpatialReference()

                    # Get the proj4 string using the WKT
                    crs.ImportFromWkt(proj)
                    proj4 = crs.ExportToProj4()

                    # Use one tif (one array) for spatial attributes
                    array = data.ReadAsArray()
                    ny, nx = array.shape
                    xs = np.arange(nx) * geom[1] + geom[0]
                    ys = np.arange(ny) * geom[5] + geom[3]

                    # Today's date for attributes
                    todays_date = dt.datetime.today()
                    today = np.datetime64(todays_date)

                    # Create Dataset
                    nco = Dataset(nc_file_name, mode="w", format="NETCDF4", clobber=True)

                    # Dimensions
                    nco.createDimension("y", ny)
                    nco.createDimension("x", nx)
                    nco.createDimension("time", None)

                    # Variables
                    y = nco.createVariable("y", np.float64, ("y",))
                    x = nco.createVariable("x", np.float64, ("x",))
                    times = nco.createVariable("time", np.int64, ("time",))
                    variable = nco.createVariable("value", np.int16,
                                                  ("time", "y", "x"),
                                                  fill_value=-9999, zlib=True)
                    variable.standard_name = "day"
                    variable.long_name = "Burn Days"

                    # Appending the CRS information
                    # Check "https://cf-trac.llnl.gov/trac/ticket/77"
                    crs = nco.createVariable("crs", "c")
                    variable.setncattr("grid_mapping", "crs")
                    crs.spatial_ref = proj4
                    crs.proj4 = proj4
                    crs.geo_transform = geom
                    crs.grid_mapping_name = "sinusoidal"
                    crs.false_easting = 0.0
                    crs.false_northing = 0.0
                    crs.longitude_of_central_meridian = 0.0
                    crs.longitude_of_prime_meridian = 0.0
                    crs.semi_major_axis = 6371007.181
                    crs.inverse_flattening = 0.0

                    # Coordinate attributes
                    x.standard_name = "projection_x_coordinate"
                    x.long_name = "x coordinate of projection"
                    x.units = "m"
                    y.standard_name = "projection_y_coordinate"
                    y.long_name = "y coordinate of projection"
                    y.units = "m"

                    # Other attributes
                    nco.title = "Burn Days"
                    nco.subtitle = "Burn Days Detection by MODIS since 1970."
                    nco.description = "The day that a fire is detected."
                    nco.date = pd.to_datetime(str(today)).strftime("%Y-%m-%d")
                    nco.projection = "MODIS Sinusoidal"
                    nco.Conventions = "CF-1.6"

                    # Variable Attrs
                    times.units = "days since 1970-01-01"
                    times.standard_name = "time"
                    times.calendar = "gregorian"
                    days = np.array([
                        self._convert_julian_date(d[0], d[1]) for d in [self._extract_date_parts(f) for f in files]
                    ])

                    # Write dimension data
                    x[:] = xs
                    y[:] = ys
                    times[:] = days

                    # One file a time, write the arrays
                    for tile_index, f in tqdm(enumerate(files), position=0, file=sys.stdout):
                        match = re.match(self._hdf_regex, f)
                        if match is None:
                            continue

                        regex_group_dict = match.groupdict()

                        ds = gdal.Open(f).GetSubDatasets()[0][0]
                        hdf = gdal.Open(ds)
                        data = hdf.GetRasterBand(1)
                        array = data.ReadAsArray()
                        year = int(regex_group_dict['year'])
                        array = self._convert_dates(array, year)
                        if array.shape == (ny, nx):
                            variable[tile_index, :, :] = array
                        else:
                            print(f + ": failed, had wrong dimensions, inserting a blank array in its place.")
                            variable[tile_index, :, :] = np.zeros((ny, nx))

                    # Done
                    nco.close()

            except Exception as e:
                # Log the error and move on to the next tile
                logging.error(f"Error processing tile {tile_id}: {str(e)}")
