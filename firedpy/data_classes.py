import datetime as dt
import os
import shutil
import sys
from datetime import datetime
from glob import glob
from typing import List, Union

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


class Base:
    """
    """
    MODIS_CRS = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    MODIS_SINUSOIDAL_PATH = os.path.join(PROJECT_DIR, "ref", "modis_grid.gpkg")
    CONUS_SHAPEFILE_PATH = os.path.join(PROJECT_DIR, 'ref', 'boundaries', 'conus.gpkg')
    DEFAULT_TILES = ["h08v04", "h09v04", "h10v04", "h11v04", "h12v04",
                     "h13v04", "h08v05", "h09v05", "h10v05", "h11v05",
                     "h12v05", "h08v06", "h09v06", "h10v06", "h11v06"]

    def __init__(self, out_dir: str, start_year: Union[None, int], end_year: Union[None, int], username: str,
                 password: str, tiles: List[str] = None):
        self._out_dir = out_dir
        self._start_year = start_year
        self._end_year = end_year
        self._username = username
        self._password = password
        self._date = datetime.today().strftime('%m-%d-%Y')
        self._cpus = os.cpu_count()
        self._land_cover_save_dir = os.path.join(out_dir, 'rasters', 'landcover')
        self._land_cover_file_root = 'lc_mosaic_'
        self._nc_save_dir = os.path.join(out_dir, 'rasters', 'burn_area', 'netcdfs')
        self._hdf_save_dir = os.path.join(out_dir, 'rasters', 'burn_area', 'hdfs')
        # self._tiles = self.DEFAULT_TILES if tiles is None else tiles

        # Initialize output directory folders and files
        self._create_paths()
        self._get_shape_files()
        print("Project Folder: " + out_dir)

    def _create_paths(self):
        sub_folders = [
            os.path.join('rasters', 'burn_area'),
            self._nc_save_dir,
            self._hdf_save_dir,
            os.path.join('rasters', 'ecoregion'),
            self._land_cover_save_dir,
            os.path.join('rasters', 'landcover', 'mosaics'),
            os.path.join('shapefiles', 'ecoregion'),
            os.path.join('tables')
        ]

        for f in [os.path.join(self._out_dir, sf) for sf in sub_folders]:
            os.makedirs(f, exist_ok=True)

    def _get_shape_files(self):
        """
        Just to grab some basic shapefiles needed for calculating statistics.
        """
        shape_dir = os.path.join(self._out_dir, "shapefiles")
        os.makedirs(shape_dir, exist_ok=True)

        files_to_copy = {
            'modis_sinusoidal_grid_world.shp': self.MODIS_SINUSOIDAL_PATH,
            'conus.shp': self.CONUS_SHAPEFILE_PATH
        }

        for dest_name, source_path in files_to_copy.items():
            dest_path = os.path.join(shape_dir, dest_name)
            if not os.path.exists(dest_path):
                shutil.copy(source_path, dest_path)

        run_conus_modis_path = os.path.join(shape_dir, 'conus_modis.shp')
        if not os.path.exists(run_conus_modis_path):
            print("Reprojecting state shapefile to MODIS Sinusoidal...")
            conus_path = os.path.join(shape_dir, 'conus.shp')
            conus = gpd.read_file(conus_path)
            modis_conus = conus.to_crs(self.MODIS_CRS)
            modis_conus.to_file(run_conus_modis_path)


class BurnData(Base):
    def __init__(self, out_dir: str, start_year: Union[None, int], end_year: Union[None, int], username: str,
                 password: str, tiles: List[str] = None):
        super().__init__(out_dir, start_year, end_year, username, password, tiles)
        self._base_sftp_folder = os.path.join('data', 'MODIS', 'C61', 'MCD64A1', 'HDF')
        self._modis_template_path = os.path.join(out_dir, 'rasters', 'mosaic_template.tif')

    @staticmethod
    # TODO: Actually look at the files in the sftp dir and build the file names from there
    def _get_tile_range(start_year, end_year):
        if start_year is None or end_year is None:
            return []
        return ["MCD64A1.A" + str(yr) for yr in range(start_year, end_year + 1)]

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

    @staticmethod
    def _convert_julien_date(julien_day, year):
        base = dt.datetime(1970, 1, 1)
        date = dt.datetime(year, 1, 1) + dt.timedelta(int(julien_day))
        days = date - base
        return days.days

    def _convert_dates(self, array, year):
        """Convert every day in an array to days since Jan 1 1970"""
        # Loop through each position with data and convert
        locs = np.where(array > 0)
        ys = locs[0]
        xs = locs[1]
        locs = [[ys[i], xs[i]] for i in range(len(xs))]
        for loc in locs:
            y = loc[0]
            x = loc[1]
            array[y, x] = self._convert_julien_date(array[y, x], year)

        return array

    def _download_files(self, sftp_client, tile: str, hdfs: List[str]) -> None:
        for hdf_file in tqdm(hdfs):
            #try:
            remote_path = self._generate_remote_hdf_path(tile, hdf_file)
            local_path = self._generate_local_hdf_path(tile, hdf_file)
            os.makedirs(os.path.dirname(local_path), exist_ok=True)

            sftp_client.get(remote_path, local_path)
            # TODO: Need more specific exception handling here
            # except Exception as e:
            #    print(e)

    def _generate_local_hdf_path(self, tile: str, hdf_name: str) -> str:
        return os.path.join(self._base_sftp_folder, tile, hdf_name)

    def _generate_remote_hdf_path(self, tile: str, hdf_name: str) -> str:
        return os.path.join(self._base_sftp_folder, tile, hdf_name)

    def _generate_local_nc_path(self, tile: str) -> str:
        return os.path.join(self._out_dir, 'rasters', 'burn_area', 'netcdfs', f"{tile}.nc")

    def _check_and_download_missing_files(self, sftp_client, tile: str, hdfs: List[str]) -> None:
        missing_files = []
        for hdf_file in hdfs:
            local_hdf_path = self._generate_local_hdf_path(tile, hdf_file)
            remote_hdf_path = self._generate_remote_hdf_path(tile, hdf_file)
            if not os.path.exists(local_hdf_path):
                missing_files.append(remote_hdf_path)
            elif not self._verify_hdf_file(local_hdf_path):
                os.remove(local_hdf_path)
                missing_files.append(remote_hdf_path)

        if not missing_files:
            return

        print(f"Missed Files: {missing_files}")
        print("trying again...")
        self._download_files(sftp_client, tile, missing_files)

    def get_burns(self, tiles: List[str], start_data: datetime, end_date: datetime):
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

        # Open the connection to the SFTP
        sftp_client = ssh_client.open_sftp()

        for tile in tiles:
            nc_file = self._generate_local_nc_path(tile)

            if os.path.exists(nc_file):
                continue

            print(f"Downloading/Checking HDF files for: {tile}")
            sftp_client.chdir(os.path.join(self._base_sftp_folder, tile))
            hdfs = [h for h in sftp_client.listdir() if ".hdf" in h]
            if not hdfs:
                print(f"No MCD64A1 Product for tile: {tile}, skipping...")

            self._download_files(sftp_client, tile, hdfs)
            self._check_and_download_missing_files(sftp_client, tile, hdfs)

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
        folders = glob(os.path.join(self._hdf_save_dir, "*"))
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

    def _write_ncs(self, tiles: List[str]):
        """
        Take in a time series of files for the MODIS burn detection dataset and
        create a singular netcdf file.
        """
        # Build the net cdfs here
        for tile_id in tiles:
            files = glob(os.path.join(self._hdf_save_dir, tile_id, "*hdf")).sort()
            # try:
            # Set file names
            # TODO: Definitely need a regex or file class here
            names = [os.path.split(f)[-1] for f in files]
            names = [f.split(".")[2] + "_" + f.split(".")[1][1:] for f in names]
            tile_id = names[0].split("_")[0]
            file_name = os.path.join(self._nc_save_dir, tile_id + ".nc")

            # Skip if it exists already
            if os.path.exists(file_name):
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

                # Get the proj4 string usign the WKT
                crs.ImportFromWkt(proj)
                proj4 = crs.ExportToProj4()

                # Use one tif (one array) for spatial attributes
                array = data.ReadAsArray()
                ny, nx = array.shape
                xs = np.arange(nx) * geom[1] + geom[0]
                ys = np.arange(ny) * geom[5] + geom[3]

                # Todays date for attributes
                todays_date = dt.datetime.today()
                today = np.datetime64(todays_date)

                # Create Dataset
                nco = Dataset(file_name, mode="w", format="NETCDF4", clobber=True)

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
                datestrings = [f[-7:] for f in names]
                dates = []
                for d in datestrings:
                    year = dt.datetime(year=int(d[:4]), month=1, day=1)
                    date = year + dt.timedelta(int(d[4:]))
                    dates.append(date)
                deltas = [d - dt.datetime(1970, 1, 1) for d in dates]
                days = np.array([d.days for d in deltas])

                # Write dimension data
                x[:] = xs
                y[:] = ys
                times[:] = days

                # One file a time, write the arrays
                tidx = 0
                for f in tqdm(files, position=0, file=sys.stdout):
                    ds = gdal.Open(f).GetSubDatasets()[0][0]
                    hdf = gdal.Open(ds)
                    data = hdf.GetRasterBand(1)
                    array = data.ReadAsArray()
                    year = int(f[-36: -32])
                    array = self._convert_dates(array, year)
                    #try:
                    variable[tidx, :, :] = array
                    #except Exception:  # TODO: Need more explicit exception handling here
                    #print(f + ": failed, probably had wrong dimensions, " + "inserting a blank array in its place.")
                    #blank = np.zeros((ny, nx))
                    #variable[tidx, :, :] = blank
                    tidx += 1

                # Done
                nco.close()

            # except Exception as e:  # TODO: More explicit exception handling
            #     file_name = os.path.join(self._nc_save_dir, tile_id + ".nc")
            #     print("Error on tile " + tile_id + ": " + str(e))
            #     print("Removing " + file_name + " and moving on.")
            #     os.remove(file_name)

