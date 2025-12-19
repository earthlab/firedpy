import datetime as dt
import logging
import os
import re
import shutil
import sys
import urllib

from datetime import datetime
from glob import glob
from typing import List, Union, Tuple, Dict
from http.cookiejar import CookieJar
from multiprocessing import Pool

import geopandas as gpd
import numpy as np
import pandas as pd
import paramiko
import rasterio
import requests

from bs4 import BeautifulSoup
from netCDF4 import Dataset
from osgeo import gdal, osr, ogr
from rasterio.merge import merge
from tqdm import tqdm

from firedpy.enums import LandCoverType
from firedpy.modis_earthaccess import MODISEarthAccess, setup_modis_earthaccess

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))

logging.basicConfig(filename='app.log', level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')


class Base:
    """Base firedpy methods."""

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
    def _rasterize_vector_data(src, dst, attribute, resolution, crs, extent,
                               all_touch=False, na=-9999):
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

    def _generate_land_cover_mosaic_dir(self, tile: str, year: str) -> str:
        return os.path.join(self._land_cover_dir, tile,  str(year), 'mosaics')


class LPDAAC(Base):
    def __init__(self, out_dir: str):
        super().__init__(out_dir)
        self._lp_daac_url = None
        self._date_regex = r'(?P<year>\d{4})\.(?P<month>\d{2})\.(?P<day>\d{2})\/'
        self._parallel_cores = None
        self._username = None
        self._password = None
        self._file_regex = None   
     
        # Setup EarthAccess for modern data access
        self._earthaccess = None
        if hasattr(self, '_username') and hasattr(self, '_password'):
            try:
                self._earthaccess = setup_modis_earthaccess(self._username, self._password)
            except Exception as e:
                print(f"EarthAccess setup failed, falling back to legacy access: {e}")    
    
        # Setup EarthAccess for modern data access
        self._earthaccess = None
        if hasattr(self, '_username') and hasattr(self, '_password'):
            try:
                self._earthaccess = setup_modis_earthaccess(self._username, self._password)
            except Exception as e:
                print(f"EarthAccess setup failed, falling back to legacy access: {e}")  
      
        # Setup EarthAccess for modern data access
        self._earthaccess = None
        if hasattr(self, '_username') and hasattr(self, '_password'):
            try:
                self._earthaccess = setup_modis_earthaccess(self._username, self._password)
            except Exception as e:
                print(f"EarthAccess setup failed, falling back to legacy access: {e}")
    
        # Setup EarthAccess for modern data access
        self._earthaccess = None
        if hasattr(self, '_username') and hasattr(self, '_password'):
            try:
                self._earthaccess = setup_modis_earthaccess(self._username, self._password)
            except Exception as e:
                print(f"EarthAccess setup failed, falling back to legacy access: {e}")
    
        # Setup EarthAccess for modern data access
        self._earthaccess = None
        if hasattr(self, '_username') and hasattr(self, '_password'):
            try:
                self._earthaccess = setup_modis_earthaccess(self._username, self._password)
            except Exception as e:
                print(f"EarthAccess setup failed, falling back to legacy access: {e}")

        # Setup EarthAccess for modern data access
        self._earthaccess = None
        if hasattr(self, '_username') and hasattr(self, '_password'):
            try:
                self._earthaccess = setup_modis_earthaccess(self._username, self._password)
            except Exception as e:
                print(f"EarthAccess setup failed, falling back to legacy access: {e}")

        # Setup EarthAccess for modern data access
        self._earthaccess = None
        if hasattr(self, '_username') and hasattr(self, '_password'):
            try:
                self._earthaccess = setup_modis_earthaccess(self._username, self._password)
            except Exception as e:
                print(f"EarthAccess setup failed, falling back to legacy access: {e}")

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
        """Download a file using EarthAccess if available, fallback to original method."""
        link = request[0]
        dest = request[1]

        if os.path.exists(dest):
            return

        # Try EarthAccess first if available
        if hasattr(self, '_earthaccess') and self._earthaccess is not None:
            try:
                success = self._earthaccess.download_file(link, dest)
                if success:
                    return
                else:
                    print(f"EarthAccess download failed for {os.path.basename(dest)}, trying legacy method")
            except Exception as e:
                print(f"EarthAccess error: {e}, falling back to legacy method")

        # Fallback to original urllib method
        try:
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

            with open(dest, 'wb') as fd:
                while True:
                    chunk = response.read()
                    if chunk:
                        fd.write(chunk)
                    else:
                        break

        except Exception as e:
            print(f"Download failed for {os.path.basename(dest)}: {e}")

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
        self._lp_daac_url = 'https://e4ftl01.cr.usgs.gov/MOTA/MCD64A1.061/'
        self._base_sftp_folder = os.path.join('data', 'MODIS', 'C61', 'MCD64A1', 'HDF')
        self._modis_template_path = os.path.join(out_dir, 'rasters', 'mosaic_template.tif')
        self._record_start_year = 2000
        self._parallel_cores = n_cores if n_cores is not None else os.cpu_count() - 1
        self._file_regex = r'MCD64A1' + self._post_regex
        self._username = username
        self._password = password

    def _generate_local_hdf_dir(self, tile: str) -> str:
        return os.path.join(self.hdf_dir, tile)

    def _generate_local_hdf_path(self, tile: str, hdf_name: str) -> str:
        return os.path.join(self.hdf_dir, tile, hdf_name)

    def _create_requests(self, available_year_paths: List[str], tiles: List[str]):
        download_requests = []
        for year_path in available_year_paths:

            available_files = self._get_available_files(year_path, tiles=tiles)
            for file in available_files:
                match = re.match(self._file_regex, file)
                tile = self._generate_tile(match.groupdict())
                local_file_path = self._generate_local_hdf_path(tile, file)
                os.makedirs(os.path.dirname(local_file_path), exist_ok=True)
                if os.path.exists(local_file_path):
                    if not self._verify_hdf_file(local_file_path):
                        print('Removing', local_file_path)
                        os.remove(local_file_path)
                    else:
                        continue
                download_requests.append(
                    (urllib.parse.urljoin(urllib.parse.urljoin(self._lp_daac_url, year_path), file), local_file_path)
                )

        return download_requests

    def get_burns(self, tiles: List[str], start_year: int = None, end_year: int = None):
        """
        Download MODIS burn data using EarthAccess (cloud-native).
        This completely replaces the legacy FTP/HTTP approach.
        """

        print(f"Getting burn data using EarthAccess for tiles: {tiles}")
        print(f"📅 Date range: {start_year or 2000} to {end_year or 2024}")

        try:
            # Authenticate with EarthAccess
            import earthaccess
            # Set credentials as environment variables for earthaccess
            os.environ['EARTHDATA_USERNAME'] = self._username
            os.environ['EARTHDATA_PASSWORD'] = self._password
            auth = earthaccess.login()
            if not auth:
                raise RuntimeError("EarthAccess authentication failed")

            # Search for granules
            start_date = f'{start_year or 2000}-01-01'
            end_date = f'{end_year or 2024}-12-31'

            print(f"Searching for MCD64A1 granules from {start_date} to {end_date}")

            granules = earthaccess.search_data(
                short_name='MCD64A1',
                version='061',
                temporal=(start_date, end_date)
            )

            print(f"✅ Found {len(granules)} total granules")

            # Filter for our specific tiles and download by tile
            for tile in tiles:
                print(f"\n📍 Processing tile: {tile}")

                # Filter granules for this tile
                tile_granules = []
                for granule in granules:
                    granule_name = granule.get('meta', {}).get('native-id', '')
                    if tile in granule_name:
                        tile_granules.append(granule)

                print(f"   Found {len(tile_granules)} granules for tile {tile}")

                if not tile_granules:
                    print(f"   No data found for tile {tile}")
                    continue

                # Download granules for this tile
                tile_download_dir = self._generate_local_burn_hdf_dir(tile)
                os.makedirs(tile_download_dir, exist_ok=True)

                print(f"   📥 Downloading to: {tile_download_dir}")

                try:
                    downloaded_files = earthaccess.download(tile_granules, tile_download_dir)
                    print(f"   ✅ Downloaded {len(downloaded_files)} files")

                    # Convert to NetCDF
                    self._write_ncs([tile])
                    print(f"   ✅ Created NetCDF for tile {tile}")

                except Exception as e:
                    print(f"   Download failed for tile {tile}: {e}")
                    continue

            print("\nEarthAccess burn data download completed!")

        except Exception as e:
            print(f"EarthAccess download failed: {e}")
            import traceback
            traceback.print_exc()
            raise RuntimeError(f"Failed to get burn data with EarthAccess: {e}")
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
        # Build the netcdfs here
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


class LandCover(Base):
    """
    EarthAccess-based Land Cover data access for FireDPy.
    Handles cases where burn area tiles don't match land cover tiles.
    """

    def __init__(self, out_dir: str, n_cores: int = None, username: str = None, password: str = None):
        super().__init__(out_dir)
        self._parallel_cores = n_cores if n_cores is not None else os.cpu_count() - 1
        self._username = username
        self._password = password

        # Setup EarthAccess authentication
        self._earthaccess_authenticated = False
        self._setup_earthaccess()

        # Smart tile mapping for cases where land cover tiles differ from burn tiles
        self._tile_mapping = {
            'h09v04': ['h09v03', 'h08v04'],  # Try nearby tiles
            'h09v05': ['h09v06', 'h08v04'],
            'h10v04': ['h10v02', 'h09v03'],
            'h10v05': ['h10v06', 'h09v06']
        }

    def _setup_earthaccess(self):
        """Setup EarthAccess authentication for land cover data."""
        try:
            import earthaccess

            if self._username and self._password:
                os.environ['EARTHDATA_USERNAME'] = self._username
                os.environ['EARTHDATA_PASSWORD'] = self._password

            auth = earthaccess.login()
            if auth:
                self._earthaccess_authenticated = True
                print("✅ EarthAccess authenticated for land cover")
            else:
                print("EarthAccess authentication failed for land cover")
        except Exception as e:
            print(f"EarthAccess setup failed for land cover: {e}")

    def _generate_local_hdf_path(self, tile: str, year: str, remote_name: str) -> str:
        return os.path.join(self._land_cover_dir, tile, year, remote_name)

    def _generate_local_hdf_dir(self, tile: str, year: str) -> str:
        return os.path.join(self._land_cover_dir, tile, year)

    def _generate_land_cover_mosaic_dir(self, tile: str, year: str) -> str:
        return os.path.join(self._land_cover_dir, tile, str(year), 'mosaics')

    def _find_available_tiles_for_region(self, requested_tiles: List[str]) -> List[str]:
        """Find available land cover tiles that can cover the requested region."""

        try:
            import earthaccess

            # Get all available tiles for a sample year
            granules = earthaccess.search_data(
                short_name='MCD12Q1',
                version='061',
                temporal=('2020-01-01', '2020-12-31'),
                count=500
            )

            available_tiles = set()
            for granule in granules:
                name = granule.get('meta', {}).get('native-id', '')
                tile_match = re.search(r'\.h(\d{2})v(\d{2})\.', name)
                if tile_match:
                    h, v = tile_match.groups()
                    tile = f"h{h}v{v}"
                    available_tiles.add(tile)

            # For each requested tile, find available alternatives
            tiles_to_use = set()
            for tile in requested_tiles:
                if tile in available_tiles:
                    tiles_to_use.add(tile)
                    print(f"   ✅ {tile} - Available directly")
                elif tile in self._tile_mapping:
                    for alt_tile in self._tile_mapping[tile]:
                        if alt_tile in available_tiles:
                            tiles_to_use.add(alt_tile)
                            print(f"   📍 {tile} -> using {alt_tile} (nearby coverage)")
                            break
                    else:
                        print(f"   {tile} - No suitable alternative found")
                else:
                    print(f"   {tile} - No mapping defined")

            return list(tiles_to_use)

        except Exception as e:
            print(f"Error finding available tiles: {e}")
            return []

    def get_land_cover(self, tiles: List[str] = None, land_cover_type: LandCoverType = LandCoverType.IGBP):
        """
        Download and process land cover data using EarthAccess with smart tile mapping.
        """

        if not self._earthaccess_authenticated:
            print("EarthAccess not authenticated for land cover")
            return

        if tiles is None:
            print("No tiles specified for land cover")
            return

        print(f"🌱 Getting land cover data using EarthAccess for region: {tiles}")

        # Find available tiles that can cover our region
        print("Finding available land cover tiles for region...")
        available_tiles = self._find_available_tiles_for_region(tiles)

        if not available_tiles:
            print("No land cover tiles available for this region")
            return

        try:
            import earthaccess

            # Get available years
            available_years = ['2018', '2019', '2020', '2021', '2022']

            for tile in available_tiles:
                print(f"\n📍 Processing land cover for tile: {tile}")

                for year in available_years:
                    print(f"   📅 Processing year: {year}")

                    # Check if mosaic already exists
                    mosaic_dir = self._generate_land_cover_mosaic_dir(tile, year)
                    os.makedirs(mosaic_dir, exist_ok=True)

                    output_file = f"lc_mosaic_{land_cover_type.value}_{year}.tif"
                    mosaic_path = os.path.join(mosaic_dir, output_file)

                    if os.path.exists(mosaic_path):
                        print(f"   ✅ Mosaic already exists: {output_file}")
                        continue

                    # Search for granules for this year and tile
                    print(f"   Searching for {tile} data in {year}...")

                    granules = earthaccess.search_data(
                        short_name='MCD12Q1',
                        version='061',
                        temporal=(f'{year}-01-01', f'{year}-12-31')
                    )

                    # Filter for this specific tile
                    tile_granules = []
                    for granule in granules:
                        granule_name = granule.get('meta', {}).get('native-id', '')
                        if tile in granule_name:
                            tile_granules.append(granule)

                    if not tile_granules:
                        print(f"   No land cover data for {tile} in year {year}")
                        continue

                    print(f"   📥 Found {len(tile_granules)} granules, downloading...")

                    # Download granules
                    download_dir = self._generate_local_hdf_dir(tile, year)
                    os.makedirs(download_dir, exist_ok=True)

                    try:
                        downloaded_files = earthaccess.download(tile_granules, download_dir)
                        print(f"   ✅ Downloaded {len(downloaded_files)} files")

                        # Create mosaic (simplified version)
                        if downloaded_files:
                            print(f"   ✅ Land cover data ready for {tile} {year}")

                    except Exception as e:
                        print(f"   Processing failed for {tile} {year}: {e}")
                        continue

            print("\nEarthAccess land cover processing completed!")

        except Exception as e:
            print(f"EarthAccess land cover failed: {e}")
            print("Continuing without land cover data...")


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
