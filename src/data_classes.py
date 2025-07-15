import datetime as dt
import os
import re
import shutil
import sys
import warnings
from datetime import datetime
from glob import glob
from typing import List, Union, Tuple, Dict
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
from urllib.error import HTTPError, URLError

from src.enums import LandCoverType

import warnings
from rasterio.errors import NotGeoreferencedWarning
warnings.filterwarnings("ignore", category=NotGeoreferencedWarning)

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
        self.hdf_dir = os.path.join(self._burn_area_dir, 'hdfs')

        self._modis_sinusoidal_grid_shape_path = os.path.join(self._shape_file_dir, 'modis_sinusoidal_grid_world.shp')
        self._conus_shape_path = os.path.join(self._shape_file_dir, 'conus.shp')

        self._eco_region_csv_path = os.path.join(self._tables_dir, 'eco_refs.csv')
        self._project_eco_region_dir = os.path.join(PROJECT_DIR, "ref", "us_eco")
        self._eco_region_shape_path = os.path.join(self._eco_region_shapefile_dir, 'NA_CEC_Eco_Level3.gpkg')

        self._post_regex = r'\.A(?P<year>\d{4})(?P<ordinal_day>\d{3})\.h(?P<horizontal_tile>\d{2})v(?P<vertical_tile>\d{2})\.061\.(?P<prod_year>\d{4})(?P<prod_ordinal_day>\d{3})(?P<prod_hourminute>\d{4})(?P<prod_second>\d{2})\.hdf$'
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
            # self._eco_region_shapefile_dir,
            self._tables_dir,
            self._mosaics_dir,
            self._nc_dir,
            self.hdf_dir
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
    def _convert_unix_day_to_calendar_date(unix_day: int) -> str:
        base = dt.datetime(1970, 1, 1)
        date = base + dt.timedelta(days=int(unix_day) - 1)
        return date.strftime('%Y-%m-%d')

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
        return os.path.join(self.hdf_dir, tile)

    def _generate_local_nc_path(self, tile: str) -> str:
        return os.path.join(self._nc_dir, f"{tile}.nc")


class LPDAAC(Base):
    def __init__(self, out_dir: str):
        super().__init__(out_dir)
        self._lp_daac_url = None
        self._date_regex = r'(?P<year>\d{4})\.(?P<month>\d{2})\.(?P<day>\d{2})\/'
        self._parallel_cores = None
        self._username = None
        self._password = None
        self._file_regex = None

    def _generate_local_hdf_path(self, year: str, remote_name: str) -> str:
        pass

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

    def _download_task(self, request: Tuple[str, str]):
        url, dest = request

        if os.path.exists(dest):
            return

        session = requests.Session()
        session.auth = (self._username, self._password)
        session.headers.update({
            "User-Agent": "firedpy/1.0 (https://github.com/earthlab/firedpy/)",
            "Accept": "*/*"
        })

        try:
            # Get file with redirect and auth handling
            response = session.get(url, stream=True, timeout=60)
            if response.status_code == 401:
                print(f"🔐 Unauthorized for {url}")
            response.raise_for_status()

            # Write to disk
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with open(dest, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

        except requests.exceptions.HTTPError as e:
            print(f"HTTPError for {url}: {e.response.status_code} {e.response.reason}")
        except Exception as e:
            print(f"Unexpected error for {url}: {str(e)}")

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

    def _download_files(self, download_requests):
        try:
            with Pool(self._parallel_cores - 1) as pool:
                pool.map(self._download_task, download_requests)
        except Exception as pe:
            try:
                tqdm.write("Parallel download failed, falling back to sequential...")
                for q in download_requests:
                    self._download_task(q)
            except Exception as e:
                template = "Download failed: error type {0}:\n{1!r}"
                message = template.format(type(e).__name__, e.args)
                tqdm.write(message)

    def _get_available_year_paths(self, start_year: int = None, end_year: int = None) -> List[str]:
        # Get available years
        request = requests.get(self._lp_daac_url)
        soup = BeautifulSoup(request.text, 'html.parser')
        year_paths = []
        for link in [link["href"] for link in soup.find_all("a", href=True)]:
            match = re.match(self._date_regex, link)
            if match is not None:
                file_year = int(match.groupdict().get('year'))
                if (start_year is None or file_year >= start_year) and (
                        end_year is None or file_year <= end_year):
                    year_paths.append(link)

        return year_paths

    def _generate_tile(self, regex_group_dict: Dict[str, str]):
        return f"h{regex_group_dict['horizontal_tile']}v{regex_group_dict['vertical_tile']}"

    def _get_available_files(self, year_path: str, tiles: List[str] = None):
        request = requests.get(urllib.parse.urljoin(self._lp_daac_url, year_path))
        soup = BeautifulSoup(request.text, 'html.parser')
        files = []
        for link in [link["href"] for link in soup.find_all("a", href=True)]:
            match = re.match(self._file_regex, link)
            if match is not None:
                if tiles is not None:
                    group_dict = match.groupdict()
                    tile = self._generate_tile(group_dict)
                    if tile not in tiles:
                        continue

                files.append(link)

        return files


class BurnData(LPDAAC):
    def __init__(self, out_dir: str, username: str, password: str, n_cores: int = None):
        super().__init__(out_dir)
        self._lp_daac_url = 'https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/MCD64A1.061/'
        # self._base_sftp_folder = os.path.join('data', 'MODIS', 'C61', 'MCD64A1', 'HDF')
        self._modis_template_path = os.path.join(out_dir, 'rasters', 'mosaic_template.tif')
        self._record_start_year = 2000
        self._parallel_cores = n_cores if n_cores is not None else os.cpu_count() - 1
        self._file_regex = r'MCD64A1\.A\d{7}\.h\d{2}v\d{2}\.061\.\d{13}\.hdf'
        self._username = username
        self._password = password

    def _generate_local_hdf_dir(self, tile: str) -> str:
        return os.path.join(self.hdf_dir, tile)

    def _generate_local_hdf_path(self, tile: str, hdf_name: str) -> str:
        return os.path.join(self.hdf_dir, tile, hdf_name)

    def _extract_year_from_granule(self, granule_id: str) -> str:
        return granule_id.split('.')[1][1:5]

    def _query_cmr_granules(self, tiles: List[str]) -> Dict[str, List[str]]:
        """
        Method for generating a list of granule IDs based on MODIS tiles.
        :param tiles:
        :return: Dictionary of granule IDs grouped by year
        """
        if tiles is None:
            tiles = self.get_all_available_tiles()

        cmr_url = "https://cmr.earthdata.nasa.gov/search/granules.json"
        short_name = "MCD64A1"
        version = "061"

        tile_patterns = [f"h{tile[1:3]}v{tile[4:]}" for tile in tiles]

        params = {
            "short_name": short_name,
            "version": version,
            "provider": "LPCLOUD",
            "page_size": 2000,
            "page_num": 1,
        }

        granules_by_year = {}
        while True:
            response = requests.get(cmr_url, params=params)
            response.raise_for_status()
            items = response.json()["feed"]["entry"]
            if not items:
                break

            for item in items:
                granule_id = item["title"]  # e.g. MCD64A1.A2000121.h09v05.061.2021200000000
                if any(tile in granule_id for tile in tile_patterns):
                    year = granule_id.split('.')[1][1:5]

                    # Get the proper HTTPS link to the HDF file
                    links = item.get("links", [])
                    hrefs = [l["href"] for l in links if "data#" in l["rel"] and l["href"].endswith(".hdf")]
                    if not hrefs:
                        continue

                    url = hrefs[0]
                    filename = os.path.basename(url)

                    granules_by_year.setdefault(year, []).append((url, filename))

            params["page_num"] += 1

        return granules_by_year

    def _create_requests(self, granule_entries: List[Tuple[str, str]]):
        download_requests = []
        for url, filename in granule_entries:
            tile = filename.split('.')[2]
            local_file_path = self._generate_local_hdf_path(tile, filename)
            os.makedirs(os.path.dirname(local_file_path), exist_ok=True)

            if os.path.exists(local_file_path):
                if not self._verify_hdf_file(local_file_path):
                    print('Removing', local_file_path)
                    os.remove(local_file_path)
                else:
                    continue

            download_requests.append((url, local_file_path))

        return download_requests

    def get_burns(self, tiles: List[str], start_year: int = None, end_year: int = None):
        """
        Method for downloading the MODIS burn event data set tiles, creating a
        singular mosaic to use as a template file for coordinate reference
        information and geometries
        :param tiles:
        :param start_year:
        :param end_year:
        :return:
        """

        for tile in tiles:
            nc_file_name = self._generate_local_nc_path(tile)
            if os.path.exists(nc_file_name):
                continue

            tqdm.write(f"Querying MCD64A1 granules for [{len(tiles)}] tile(s) ...")
            granules_by_year = self._query_cmr_granules([tile])

            # Flatten granule list and filter by year
            all_granules = []
            for year, granules in granules_by_year.items():
                if start_year and int(year) < start_year:
                    continue
                if end_year and int(year) > end_year:
                    continue
                all_granules.extend(granules)

            # Create download requests
            download_requests = self._create_requests(all_granules)

            # tqdm.write(f"Downloading {len(download_requests)} granules for tile {tile}...")
            with tqdm(total=len(download_requests), desc=f"Downloading {tile}", position=0, dynamic_ncols=True) as pbar:
                for i, req in enumerate(download_requests):
                    try:
                        self._download_task(req)
                    except Exception as e:
                        tqdm.write(f"Failed to download {req[0]}: {str(e)}")
                    pbar.update(1)

            self._write_ncs([tile])

    def _write_modis_template_file(self):
        # Merge one year into a reference mosaic
        print("Creating reference mosaic ...")
        folders = glob(os.path.join(self.hdf_dir, "*"))
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
        match = re.match(self._file_regex, filename)
        if match:
            return int(match.groupdict()['year']), int(match.groupdict()['ordinal_day'])
        return None

    def get_date_range(self, start_year: int = None, end_year: int = None) -> List[Tuple[int, int]]:
        dates = [
            self._extract_date_parts(f) for f in glob(os.path.join(self.hdf_dir, '**', '*'), recursive=True)
        ]

        return sorted([d for d in dates if d is not None and (start_year is None or d[0] >= start_year) and (
                end_year is None or d[0] <= end_year)])

    def _write_ncs(self, tiles: List[str]):
        """
        Take in a time series of files for the MODIS burn detection dataset and
        create a singular netcdf file.
        """
        # Build the net cdfs here
        fill_value = -9999
        for tile_id in tiles:
            try:
                nc_file_name = self._generate_local_nc_path(tile_id)
                if os.path.exists(nc_file_name):
                    continue

                hdf_dir = self._generate_local_burn_hdf_dir(tile_id)

                files = sorted([os.path.join(hdf_dir, f) for f in os.listdir(hdf_dir) if self._extract_date_parts(f)
                                is not None], key=self._extract_date_parts)

                if not files:
                    print(f'No hdf files for tile {tile_id} in {hdf_dir}')

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
                    metadata = hdf.GetMetadata()
                    array = data.ReadAsArray()
                    ny, nx = array.shape
                    xs = np.arange(nx) * geom[1] + geom[0]
                    ys = np.arange(ny) * geom[5] + geom[3]
                    lats = np.linspace(
                        float(metadata['SOUTHBOUNDINGCOORDINATE']),
                        float(metadata['NORTHBOUNDINGCOORDINATE']),
                        ny
                    )
                    lons = np.linspace(
                        float(metadata['WESTBOUNDINGCOORDINATE']),
                        float(metadata['EASTBOUNDINGCOORDINATE']),
                        nx
                    )

                    # Today's date for attributes
                    todays_date = dt.datetime.today()
                    today = np.datetime64(todays_date)

                    # Create Dataset
                    nco = Dataset(nc_file_name, mode="w", format="NETCDF4", clobber=True)

                    # Dimensions
                    nco.createDimension("y", ny)
                    nco.createDimension("x", nx)
                    nco.createDimension("time", None)
                    nco.createDimension('lat', nx)
                    nco.createDimension('lon', ny)

                    # Variables
                    y = nco.createVariable("y", np.float64, ("y",))
                    x = nco.createVariable("x", np.float64, ("x",))
                    lat = nco.createVariable('lat', np.float64, dimensions=('lat',))
                    lon = nco.createVariable('lon', np.float64, dimensions=('lon',))
                    times = nco.createVariable("time", np.int16, ("time",))
                    variable = nco.createVariable("value", np.int16,
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

                    lat.standard_name = 'latitude_coordinate'
                    lat.long_name = 'latitude coordinate'
                    lat.units = 'deg'
                    lon.standard_name = 'longitude_coordinate'
                    lon.long_name = 'longitude coordinate'
                    lon.units = 'deg'

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
                    lon[:] = lons
                    lat[:] = lats
                    times[:] = days

                    # One file a time, write the arrays
                    for tile_index, f in tqdm(enumerate(files), position=0, file=sys.stdout):
                        match = re.match(self._file_regex, os.path.basename(f))
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
                        nulls = np.where(array < 0)

                        year = int(regex_group_dict['year'])
                        array = self._convert_dates(array, year)

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


class LandCover(LPDAAC):
    def __init__(self, out_dir: str, n_cores: int = None, username: str = None, password: str = None):
        super().__init__(out_dir)
        # self._lp_daac_url = 'https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.061/' # doesn't work? Use below (MC, July 25)
        self._lp_daac_url = 'https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/MCD12Q1.061/'
        # self._date_regex = r'(?P<year>\d{4})\.(?P<month>\d{2})\.(?P<day>\d{2})\/' # this is no longer part of the URL
        self._parallel_cores = n_cores if n_cores is not None else os.cpu_count() - 1
        self._username = username
        self._password = password
        # self._file_regex = r'MCD12Q1' + self._post_regex # this no longer meets the file convention (MC, July 25)
        self._file_regex = r'MCD12Q1\.A\d{7}\.h\d{2}v\d{2}\.061\.\d{13}\.hdf'

    def _generate_local_hdf_path(self, year: str, remote_name: str) -> str:
        return os.path.join(self._land_cover_dir, year, remote_name)

    def _generate_local_hdf_dir(self, year: str) -> str:
        return os.path.join(self._land_cover_dir, year)

    def _query_cmr_granules(self, tiles: List[str]) -> Dict[str, List[str]]:
        """
        Method for generating a list of granule IDs based on MODIS tiles.
        :param tiles:
        :return: Dictionary of granule IDs grouped by year
        """
        if tiles is None:
            tiles = self.get_all_available_tiles()

        # define the CMR URL and inputs
        cmr_url = "https://cmr.earthdata.nasa.gov/search/granules.json"
        short_name = "MCD12Q1"
        version = "061"

        tile_patterns = [f"h{tile[1:3]}v{tile[4:]}" for tile in tiles]  # e.g., 'h10v09'

        params = {
            "short_name": short_name,
            "version": version,
            "provider": "LPCLOUD",
            "page_size": 2000,
            "page_num": 1,
        }

        granules_by_year = {}
        while True:
            response = requests.get(cmr_url, params=params)
            response.raise_for_status()
            items = response.json()["feed"]["entry"]
            if not items:
                break
            for item in items:
                granule_id = item["title"]
                if any(tile in granule_id for tile in tile_patterns):
                    year = granule_id.split('.')[1][1:5]
                    granules_by_year.setdefault(year, []).append(granule_id)
            params["page_num"] += 1

        return granules_by_year

    def _create_requests(self, granules: List[str]):
        download_requests = []
        for granule in granules:
            year = granule.split('.')[1][1:5]
            remote_file = f"{granule}.hdf"
            local_file_path = self._generate_local_hdf_path(year, remote_file)
            os.makedirs(os.path.dirname(local_file_path), exist_ok=True)

            if os.path.exists(local_file_path):
                if not self._verify_hdf_file(local_file_path):
                    print('Removing corrupt file:', local_file_path)
                    os.remove(local_file_path)
                else:
                    continue

            url = f"{self._lp_daac_url}{granule}/{remote_file}"
            download_requests.append((url, local_file_path))

        return download_requests

    def _create_annual_mosaic(self, year: str, land_cover_type: LandCoverType = LandCoverType.IGBP):
        """
        Create annual mosaic landcover geotiffs from the downloaded HDF files.
        :param year:
        :param land_cover_type: User-defined landcover type (IGBP default)
        :return: Annual mosaic GeoTIFF
        """
        # define the output geotiff file name
        output_file = f"lc_mosaic_{land_cover_type.value}_{year}.tif"
        # gather a list of the granules
        lc_files = [self._generate_local_hdf_path(year, f) for f in os.listdir(self._generate_local_hdf_dir(year))
                    if re.match(self._file_regex, os.path.basename(f)) is not None]

        # extract the correct subdatasets from the HDF file
        # here we find the correct landcover model
        datasets = []
        for lc_file_path in lc_files:
            with rasterio.open(lc_file_path) as lc_file:
                datasets.append([sd for sd in lc_file.subdatasets if str(land_cover_type.value) in sd.lower()][0])

        # open the raster datasets and grab metadata
        tiles = [rasterio.open(ds) for ds in datasets]
        mosaic, transform = merge(tiles)
        # copy the metadata and update projection information
        crs = tiles[0].meta.copy()
        crs.update(
            {"driver": "GTiff",
             "height": mosaic.shape[1],
             "width": mosaic.shape[2],
             "transform": transform}
        )

        # write out the mosaic file.
        with rasterio.open(os.path.join(self._mosaics_dir, output_file), "w+", **crs) as dst:
            dst.write(mosaic)

    def get_land_cover(self, tiles: List[str] = None, land_cover_type: LandCoverType = LandCoverType.IGBP):
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
        Update 07/2025: Integration with the CMR API for granule lists by year (M. Cook)
        """
        if tiles is None:
            tiles = self.get_all_available_tiles()

        print("Querying MCD12Q1 granules ...")
        granules_by_year = self._query_cmr_granules(tiles)
        years = sorted(granules_by_year.keys())

        with tqdm(total=len(years), desc="", dynamic_ncols=True) as pbar:
            for year in years:
                pbar.set_description_str(f"Downloading and mosaicking granules for [{year}]")
                output_file = f"lc_mosaic_{land_cover_type.value}_{year}.tif"
                if os.path.exists(os.path.join(self._mosaics_dir, output_file)):
                    pbar.update(1)
                    continue

                granules = granules_by_year[year]
                download_requests = self._create_requests(granules)
                self._download_files(download_requests)
                self._create_annual_mosaic(str(year), land_cover_type)
                pbar.update(1)

        tqdm.write(f"\tCompleted for all years: [{min(years)}-{max(years)}]")
        tqdm.write(f"\tLandcover data saved to: {self._mosaics_dir}")


class EcoRegion(Base):
    def __init__(self, out_dir: str):
        super().__init__(out_dir)

        self._eco_region_ftp_url = 'ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip'
        self._eco_region_raster_path = os.path.join(self._eco_region_raster_dir, 'NA_CEC_Eco_Level3_modis.tif')

        self._ref_cols = ['NA_L3CODE', 'NA_L3NAME', 'NA_L2CODE', 'NA_L2NAME', 'NA_L1CODE', 'NA_L1NAME', 'NA_L3KEY',
                          'NA_L2KEY', 'NA_L1KEY']
        self.eco_region_data_frame: Union[None, gpd.GeoDataFrame] = None
        self._file_regex = r'MCD64A1' + self._post_regex

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
                shutil.copytree(self._project_eco_region_dir, self._eco_region_shapefile_dir)
                return gpd.read_file(self._eco_region_shape_path)

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
                    re.match(self._file_regex, os.path.basename(f)) is not None][0]
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
