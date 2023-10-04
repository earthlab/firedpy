import datetime as dt
import os
import re
import shutil
import sys
from datetime import datetime
from glob import glob
from typing import List, Union, Tuple
import logging
import urllib
from http.cookiejar import CookieJar
from multiprocessing import Pool

import geopandas as gpd
import numpy as np
import pandas as pd
import paramiko
import rasterio
from netCDF4 import Dataset
from osgeo import gdal, osr, ogr
from rasterio.merge import merge
from tqdm import tqdm
import requests
from bs4 import BeautifulSoup

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
        self._eco_region_shapefile_dir = os.path.join(self._shape_file_dir, 'eco_region')
        self._tables_dir = os.path.join(out_dir, 'tables')

        self._mosaics_dir = os.path.join(self._land_cover_dir, 'mosaics')
        self._nc_dir = os.path.join(self._burn_area_dir, 'netcdfs')
        self._hdf_dir = os.path.join(self._burn_area_dir, 'hdfs')

        self._modis_sinusoidal_grid_shape_path = os.path.join(self._shape_file_dir, 'modis_sinusoidal_grid_world.shp')
        self._conus_shape_path = os.path.join(self._shape_file_dir, 'conus.shp')

        self._burn_hdf_regex = r'MCD64A1\.A(?P<year>\d{4})(?P<ordinal_day>\d{3})\.h(?P<horizontal_tile>\d{2})v(?P<vertical_tile>\d{2})\.061\.(?P<prod_year>\d{4})(?P<prod_ordinal_day>\d{3})(?P<prod_hourminute>\d{4})(?P<prod_second>\d{2})\.hdf$'

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
            self._modis_sinusoidal_grid_shape_path: self.MODIS_SINUSOIDAL_PATH,
            self._conus_shape_path: self.CONUS_SHAPEFILE_PATH
        }

        for dest_path, source_path in files_to_copy.items():
            if not os.path.exists(dest_path):
                shutil.copy(source_path, dest_path)

    @staticmethod
    def _convert_ordinal_to_unix_day(year: int, ordinal_day: int) -> int:
        base = dt.datetime(1970, 1, 1)
        date = dt.datetime(year, 1, 1) + dt.timedelta(int(ordinal_day - 1))
        return (date - base).days

    @staticmethod
    def _rasterize_vector_data(src, dst, attribute, resolution, crs, extent, all_touch=False, na=-9999):
        """Rasterizes input vector data"""
        # Open shapefile, retrieve the layer
        src_data = ogr.Open(src)
        layer = src_data.GetLayer()

        # Use transform to derive coordinates and dimensions
        xmin = extent[0]
        ymin = extent[1]
        xmax = extent[2]
        ymax = extent[3]

        # Create the target raster layer
        cols = int((xmax - xmin) / resolution)
        rows = int((ymax - ymin) / resolution) + 1
        trgt = gdal.GetDriverByName("GTiff").Create(dst, cols, rows, 1,
                                                    gdal.GDT_Float32)
        trgt.SetGeoTransform((xmin, resolution, 0, ymax, 0, -resolution))

        # Add crs
        refs = osr.SpatialReference()
        refs.ImportFromWkt(crs)
        trgt.SetProjection(refs.ExportToWkt())

        # Set no value
        band = trgt.GetRasterBand(1)
        band.SetNoDataValue(na)

        # Set options
        if all_touch is True:
            ops = ["-at", "ATTRIBUTE=" + attribute]
        else:
            ops = ["ATTRIBUTE=" + attribute]

        # Finally rasterize
        gdal.RasterizeLayer(trgt, [1], layer, options=ops)
        # Close target an source rasters
        del trgt
        del src_data

    def _convert_dates(self, array, year) -> np.array:
        """Convert every day in an array to days since Jan 1 1970"""
        # Loop through each position with data and convert
        ys, xs = np.where(array > 0)

        for y, x in zip(ys, xs):
            array[y, x] = self._convert_ordinal_to_unix_day(year, array[y, x])

        return array

    def _generate_local_burn_hdf_dir(self, tile: str) -> str:
        return os.path.join(self._hdf_dir, tile)


class BurnData(Base):
    def __init__(self, out_dir: str):
        super().__init__(out_dir)
        self._base_sftp_folder = os.path.join('data', 'MODIS', 'C61', 'MCD64A1', 'HDF')
        self._modis_template_path = os.path.join(out_dir, 'rasters', 'mosaic_template.tif')
        self._record_start_year = 2000

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
                if os.path.exists(local_path):
                    continue
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
            raise IOError(f'Error downloading burn data: max retries exceeded ({max_retries}). Files not downloaded or '
                          f'not able to open: {retries}')

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
                match = re.match(self._burn_hdf_regex, hdf_file)
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
        match = re.match(self._burn_hdf_regex, filename)
        if match:
            return int(match.groupdict()['year']), int(match.groupdict()['ordinal_day'])
        return None

    def _write_ncs(self, tiles: List[str]):
        """
        Take in a time series of files for the MODIS burn detection dataset and
        create a singular netcdf file.
        """
        # Build the net cdfs here
        fill_value = -9999
        for tile_id in tiles:
            try:
                hdf_dir = self._generate_local_burn_hdf_dir(tile_id)

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
                    variable = nco.createVariable("value", np.int32,
                                                  ("time", "y", "x"),
                                                  fill_value=fill_value, zlib=True)
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
                        self._convert_ordinal_to_unix_day(d[0], d[1]) for d in [self._extract_date_parts(f) for f in
                                                                                files]
                    ])

                    # Write dimension data
                    x[:] = xs
                    y[:] = ys
                    times[:] = days

                    # One file a time, write the arrays
                    for tile_index, f in tqdm(enumerate(files), position=0, file=sys.stdout):
                        match = re.match(self._burn_hdf_regex, os.path.basename(f))
                        if match is None:
                            continue

                        regex_group_dict = match.groupdict()

                        try:
                            ds = gdal.Open(f).GetSubDatasets()[0][0]
                            hdf = gdal.Open(ds)
                        except Exception as e:
                            print(f'Could not open {f} for building ncdf: {str(e)}')
                            variable[tile_index, :, :] = np.full((ny, nx), fill_value)
                            continue

                        data = hdf.GetRasterBand(1)
                        array = data.ReadAsArray()
                        year = int(regex_group_dict['year'])
                        array = self._convert_dates(array, year)

                        nulls = np.where(array < 0)
                        array[nulls] = fill_value

                        if array.shape == (ny, nx):
                            variable[tile_index, :, :] = array
                        else:
                            print(f + ": failed, had wrong dimensions, inserting a blank array in its place.")
                            variable[tile_index, :, :] = np.full((ny, nx), fill_value)

                    # Done
                    nco.close()

            except Exception as e:
                # Log the error and move on to the next tile
                logging.error(f"Error processing tile {tile_id}: {str(e)}")


class LandCover(Base):
    def __init__(self, out_dir: str):
        super().__init__(out_dir)
        self._lp_daac_url = 'https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.061/'
        self._date_regex = r'(?P<year>\d{4})\.(?P<month>\d{2})\.(?P<day>\d{2})\/'
        self._hdf_regex = r'MCD12Q1\.A(?P<year>\d{4})(?P<ordinal_day>\d{3})\.h(?P<horizontal_tile>\d{2})v(?P<vertical_tile>\d{2})\.061\.(?P<prod_year>\d{4})(?P<prod_ordinal_day>\d{3})(?P<prod_hourminute>\d{4})(?P<prod_second>\d{2})\.hdf$'
        self._core_count = os.cpu_count()

    @staticmethod
    def get_all_available_tiles() -> List[str]:
        with paramiko.SSHClient() as ssh_client:
            ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh_client.connect(hostname="fuoco.geog.umd.edu",
                               username="fire", password="burnt")
            print("Connected to 'fuoco.geog.umd.edu' ...")
            # Open the connection to the SFTP
            with ssh_client.open_sftp() as sftp_client:
                sftp_client.chdir('/data/MODIS/C61/MCD64A1/HDF')
                return sftp_client.listdir()

    @staticmethod
    def _is_file_valid(in_file: str) -> bool:
        try:
            rasterio.open(in_file, 'w+', driver='HDF4Image')
            return True
        except Exception as _:
            return False

    def _get_available_year_paths(self) -> List[str]:
        # Get available years
        request = requests.get(self._lp_daac_url)
        soup = BeautifulSoup(request.text, 'html.parser')
        year_paths = []
        for link in [link["href"] for link in soup.find_all("a", href=True)]:
            match = re.match(self._date_regex, link)
            if match is not None:
                year_paths.append(link)

        return year_paths

    def _get_available_files(self, year_path: str, tiles: List[str] = None):
        request = requests.get(urllib.parse.urljoin(self._lp_daac_url, year_path))
        soup = BeautifulSoup(request.text, 'html.parser')
        files = []
        for link in [link["href"] for link in soup.find_all("a", href=True)]:
            match = re.match(self._hdf_regex, link)
            if match is not None:

                if tiles is not None:
                    group_dict = match.groupdict()
                    tile = f"h{group_dict['horizontal_tile']}v{group_dict['vertical_tile']}"
                    if tile not in tiles:
                        continue

                files.append(link)

        return files

    def _generate_local_hdf_path(self, year: str, remote_name: str) -> str:
        return os.path.join(self._land_cover_dir, year, remote_name)

    def _generate_local_hdf_dir(self, year: str) -> str:
        return os.path.join(self._land_cover_dir, year)

    def _download_task(self, request: Tuple[str, str]):
        link = request[0]
        dest = request[1]

        if os.path.exists(dest):
            return

        pm = urllib.request.HTTPPasswordMgrWithDefaultRealm()
        pm.add_password(None, "https://urs.earthdata.nasa.gov", self._username, self._password)
        cookie_jar = CookieJar()
        opener = urllib.request.build_opener(
            urllib.request.HTTPBasicAuthHandler(pm),
            urllib.request.HTTPCookieProcessor(cookie_jar)
        )
        urllib.request.install_opener(opener)
        myrequest = urllib.request.Request(link)
        response = urllib.request.urlopen(myrequest)
        response.begin()
        with open(dest, 'wb') as fd:
            while True:
                chunk = response.read()
                if chunk:
                    fd.write(chunk)
                else:
                    break

    def _download_files(self, download_requests):
        try:
            with Pool(int(self._core_count / 2)) as pool:
                for _ in tqdm(pool.imap_unordered(self._download_task, download_requests),
                              total=len(download_requests)):
                    pass

        except Exception as pe:
            try:
                _ = [self._download_task(q) for q in tqdm(download_requests, position=0, file=sys.stdout)]
            except Exception as e:
                template = "Download failed: error type {0}:\n{1!r}"
                message = template.format(type(e).__name__, e.args)
                print(message)

    def _create_annual_mosaic(self, year: str, land_cover_type: int = 1):
        output_file = f"lc_mosaic_{land_cover_type}_{year}.tif"
        if os.path.exists(output_file):
            return

        # Filter available files for the requested tiles
        lc_files = [self._generate_local_hdf_path(year, f) for f in os.listdir(self._generate_local_hdf_dir(year))
                    if re.match(self._hdf_regex, f) is not None]

        # Use the sub-dataset name to get the right land cover type
        datasets = []
        for lc_file_path in lc_files:
            with rasterio.open(lc_file_path) as lc_file:
                datasets.append([sd for sd in lc_file.subdatasets if land_cover_type in sd.lower()][0])

        # Create pointers to the chosen land cover type
        tiles = [rasterio.open(ds) for ds in datasets]

        # Mosaic them together
        mosaic, transform = merge(tiles)

        # Get coordinate reference information
        crs = tiles[0].meta.copy()
        crs.update({"driver": "GTIFF",
                    "height": mosaic.shape[1],
                    "width": mosaic.shape[2],
                    "transform": transform})

        # Save mosaic file
        with rasterio.open(os.path.join(self._mosaics_dir, output_file), "w+", **crs) as dst:
            dst.write(mosaic)

    def get_land_cover(self, tiles: List[str] = None, land_cover_type: int = 1):
        """
        A method to download and process land cover data from  NASA's Land
        Processes Distributed Active Archive Center, which is an Earthdata
        thing. You"ll need register for a username and password, but that"s
        free. Fortunately, there is a tutorial on how to get this data:

        https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python

        sample citation for later:
           ASTER Mount Gariwang image from 2018 was retrieved from
           https://lpdaac.usgs.gov, maintained by the NASA EOSDIS Land
           Processes Distributed Active Archive Center (LP DAAC) at the USGS
           Earth Resources Observation and Science (EROS) Center, Sioux Falls,
           South Dakota. 2018, https://lpdaac.usgs.gov/resources/data-action/
           aster-ultimate-2018-winter-olympics-observer/.

        Update 10/2020: Python workflow updated to 3.X specific per documentation

        """
        if tiles is None:
            tiles = self.get_all_available_tiles()

        available_year_paths = self._get_available_year_paths()

        download_requests = []
        for year_path in available_year_paths:
            year = re.match(self._date_regex, year_path).groupdict()['year']

            available_files = self._get_available_files(year_path, tiles=tiles)
            for file in available_files:
                local_file_path = self._generate_local_hdf_path(year, file)
                if os.path.exists(local_file_path):
                    if not self._is_file_valid(local_file_path):
                        os.remove(local_file_path)
                    else:
                        continue
                download_requests.append(
                    (urllib.parse.urljoin(urllib.parse.urljoin(self._lp_daac_url, year_path), file), local_file_path)
                )

        self._download_files(download_requests)

        print("Mosaicking/remosaicking land cover tiles...")
        for year_path in tqdm(available_year_paths, position=0, file=sys.stdout):
            year = re.match(self._date_regex, year_path).groupdict()['year']
            self._create_annual_mosaic(year, land_cover_type)

        # Print location
        print(f"Land cover data saved to {self._mosaics_dir}")


class EcoRegion(Base):
    def __init__(self, out_dir: str):
        super().__init__(out_dir)

        self._eco_region_ftp_url = 'ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip'
        self._project_eco_region_path = os.path.join(PROJECT_DIR, "ref", "us_eco", "NA_CEC_Eco_Level3.shp")
        self._eco_region_shape_path = os.path.join(self._eco_region_shapefile_dir, 'NA_CEC_Eco_Level3.gpkg')
        self._eco_region_csv_path = os.path.join(self._tables_dir, 'eco_refs.csv')
        self._eco_region_raster_path = os.path.join(self._eco_region_raster_dir, 'NA_CEC_Eco_Level3_modis.tif')

        self._ref_cols = ['NA_L3CODE', 'NA_L3NAME', 'NA_L2CODE', 'NA_L2NAME', 'NA_L1CODE', 'NA_L1NAME', 'NA_L3KEY',
                          'NA_L2KEY', 'NA_L1KEY']
        self.eco_region_data_frame: Union[None, gpd.GeoDataFrame] = None

    @staticmethod
    def _normalize_string(string) -> str:
        """
        # Character cases are inconsistent between I,II and III,IV levels
        """

        def capitalize_special(s):
            """Handles capitalization for special characters within a string."""
            if "/" in s:
                s = "/".join([segment.capitalize() if segment.upper() != 'USA' else segment.upper() for segment in
                              s.split("/")])
            if "-" in s:
                s = "-".join([segment.title() if segment.upper() != 'USA' else segment.upper() for segment in
                              s.split("-")])
            return s

        # Split string into words
        words = string.split()

        # Create a new list with words formatted according to rules
        formatted_words = []
        for word in words:
            # Special handling for "USA"
            if word.upper() == "USA":
                formatted_words.append("USA")
            # Special handling for "and"
            elif word.lower() == "and":
                formatted_words.append("and")
            # Capitalization with special character handling for other words
            else:
                formatted_words.append(capitalize_special(word.capitalize()))

        # Join and return the formatted words as a single string
        return " ".join(formatted_words)

    def _read_eco_region_file(self) -> gpd.GeoDataFrame:
        """
       Update 02/2021: EPA FTP site is glitchy, adding local file in "ref"
       If download fails, use local file
       """
        if not os.path.exists(self._eco_region_shape_path):
            try:
                eco = gpd.read_file(self._eco_region_ftp_url)
                eco.crs = {"init": "epsg:5070"}
                eco.to_file(self._eco_region_shape_path)
                return gpd.read_file(self._eco_region_shape_path)
            except Exception as e:
                print("Failed to connect to EPA ftp site: using local file ...")
                shutil.copy(self._project_eco_region_path, self._eco_region_shape_path)
                return gpd.read_file(self._eco_region_shape_path)

    def create_eco_region_raster(self, tiles: List[str]):
        if self.eco_region_data_frame is None:
            self.get_eco_region()

        # We need something with the correct geometry
        template1 = gpd.read_file(self._modis_sinusoidal_grid_shape_path)

        # Getting the extent regardless of existing files from other runs
        template1["h"] = template1["h"].apply(lambda x: "{:02d}".format(x))
        template1["v"] = template1["v"].apply(lambda x: "{:02d}".format(x))
        template1["tile"] = "h" + template1["h"] + "v" + template1["v"]
        template1 = template1[template1["tile"].isin(tiles)]

        # We can use this to query which tiles are needed for coordinates
        bounds = template1.geometry.bounds
        minx = min(bounds["minx"])
        miny = min(bounds["miny"])
        maxx = max(bounds["maxx"])
        maxy = max(bounds["maxy"])
        minx_tile = template1["tile"][bounds["minx"] == minx].iloc[0]
        miny_tile = template1["tile"][bounds["miny"] == miny].iloc[0]
        maxx_tile = template1["tile"][bounds["maxx"] == maxx].iloc[0]
        maxy_tile = template1["tile"][bounds["maxy"] == maxy].iloc[0]
        extent_tiles = [minx_tile, miny_tile, maxx_tile, maxy_tile]

        # If these aren't present, I say just go ahead and download
        exts = []
        projection = None
        xres = None
        for tile in extent_tiles:
            burn_dir = self._generate_local_burn_hdf_dir(tile)
            if not os.path.exists(burn_dir):
                raise FileNotFoundError(f'Burn HDFs do not exist for tile {tile} at {burn_dir}. Use the BurnData class '
                                        f'to download this data before rasterizing the eco region file')

            file = [os.path.join(burn_dir, f) for f in glob(os.path.join(burn_dir, "*")) if
                    re.match(self._burn_hdf_regex, os.path.basename(f)) is not None][0]
            file_pointer = gdal.Open(file)
            dataset_pointer = file_pointer.GetSubDatasets()[0][0]
            ds = gdal.Open(dataset_pointer)
            geom = ds.GetGeoTransform()
            ulx, xres, xskew, uly, yskew, yres = geom
            lrx = ulx + (ds.RasterXSize * xres)
            lry = uly + (ds.RasterYSize * yres)
            exts.append([ulx, lry, lrx, uly])
            projection = ds.GetProjection()

        extent = [exts[0][0], exts[1][1], exts[2][2], exts[3][3]]
        attribute = "US_L3CODE"
        self._rasterize_vector_data(self.eco_region_data_frame, self._eco_region_raster_path, attribute, xres,
                                    projection, extent)

    def get_eco_region(self):
        # Omernick's Ecoregions - EPA North American Albers
        eco = self._read_eco_region_file()

        # Create a reference table for ecoregions
        eco_ref = eco[self._ref_cols].drop_duplicates()
        eco_ref = eco_ref.applymap(self._normalize_string)

        eco_ref.to_csv(self._eco_region_csv_path, index=False)

        self.eco_region_data_frame = eco_ref
