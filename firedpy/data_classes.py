import datetime as dt
import logging
import os
import re
import shutil
import sys
import traceback
import urllib

from concurrent.futures import as_completed, ThreadPoolExecutor
from datetime import datetime
from glob import glob
from http.cookiejar import CookieJar
from multiprocessing import Pool
from pathlib import Path, PosixPath
from typing import List, Tuple, Dict

import earthaccess
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

from firedpy import DATA_DIR
from firedpy.enums import EcoRegionType, LandCoverType
from firedpy.modis_earthaccess import MODISEarthAccess
from firedpy.utilities.spatial import (
    country_to_tiles,
    hdf4_to_geotiff,
    shape_to_tiles,
    tiles_to_points
)

gdal.UseExceptions()
logger = logging.getLogger(__name__)


ATTR_DESCS = {
    "ecoregions": {
        "na": "North American ecoregions (Omernick, 1987)",
        "world": "World Terrestrial Ecoregions (World Wildlife Fund)"
    },
    "landcover": {
        "1": "Annual International Geosphere-Biosphere Programme (IGBP)",
        "2": "Annual University of Maryland (UMD)",
        "3": "Annual Leaf Area Index (LAI)",
        "4": "Annual BIOME-Biogeochemical Cycles (BGC)",
        "5": "Annual Plant Functional Types (PFT)"
    }
}


class Base(MODISEarthAccess):
    """Base firedpy methods."""

    MODIS_CRS = ("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 "
                 "+b=6371007.181 +units=m +no_defs")
    MODIS_SINUSOIDAL_PATH = DATA_DIR.joinpath("modis_grid.gpkg")
    CONUS_SHAPEFILE_PATH = DATA_DIR.joinpath("boundaries", "conus.gpkg")
    DEFAULT_TILES = ["h08v04", "h09v04", "h10v04", "h11v04", "h12v04",
                     "h13v04", "h08v05", "h09v05", "h10v05", "h11v05",
                     "h12v05", "h08v06", "h09v06", "h10v06", "h11v06"]

    def __init__(self, project_directory, n_cores=0):
        """Initialize Base Firedpy object.

        Parameters
        ----------
        project_directory : str | pathlib.PosixPath
            Target firedpy file output directory.
        """
        super().__init__()

        # Set number of cpu cores
        n_cores = int(n_cores)
        if n_cores is not None and n_cores != 0:
            self.n_cores = n_cores
        else:
            self.n_cores = os.cpu_count()

        # Set up output directory paths
        project_directory = os.path.expanduser(project_directory)
        os.makedirs(project_directory, exist_ok=True)
        project_directory = Path(project_directory).expanduser()
        self.project_directory = project_directory
        self.date = datetime.today().strftime("%m-%d-%Y")
        self.raster_dir = project_directory.joinpath("rasters")
        self.shape_dir = project_directory.joinpath("shapefiles")
        self.burn_area_dir = self.raster_dir.joinpath("burn_area")
        self.land_cover_dir = self.raster_dir.joinpath("land_cover")
        self.eco_region_raster_dir = self.raster_dir.joinpath("eco_region")
        self.eco_region_shape_dir = self.shape_dir.joinpath("eco_region")
        self.tables_dir = project_directory.joinpath("tables")
        self.nc_dir = self.burn_area_dir.joinpath("netcdfs")
        self.hdf_dir = self.burn_area_dir.joinpath("hdfs")
        self._modis_sinusoidal_grid_shape_path = self.shape_dir.joinpath(
            "modis_sinusoidal_grid_world.shp"
        )
        self.conus_shape_path = self.shape_dir.joinpath("conus.shp")
        self.eco_region_csv_path = self.tables_dir.joinpath("eco_refs.csv")
        self.project_eco_region_dir = DATA_DIR.joinpath("us_eco")
        self.eco_region_shape_path = self.eco_region_shape_dir.joinpath(
            "NA_CEC_Eco_Level3.gpkg"
        )

        # Can we shorten this or use another method?
        # This is used both to ensure the file format matches and to extract
        # The year and day from the file name
        post_regex = (
            r"\.A(?P<year>\d{4})"
            r"(?P<ordinal_day>\d{3})\.h"
            r"(?P<horizontal_tile>\d{2})"
            r"v(?P<vertical_tile>\d{2})\.061\."
            r"(?P<prod_year>\d{4})"
            r"(?P<prod_ordinal_day>\d{3})"
            r"(?P<prod_hourminute>\d{4})"
            r"(?P<prod_second>\d{2})\.hdf$"
        )
        self._file_regex = self._file_regex = r"MCD64A1" + post_regex

        # Initialize output directory folders and files
        self._initialize_save_dirs()
        self._get_shape_files()

    def _authenticate(self):
        # This will use a config file or a prompt if one isn't available
        auth = earthaccess.login(strategy="all")
        if not auth:
            raise RuntimeError("EarthAccess authentication failed")
        return auth

    def _copy_cec_file(self):
        src = DATA_DIR.joinpath("ec_eco", "NA_CEC_Eco_Level3.gpkg")
        return self._copy_file(src, self.eco_region_shape_path)

    @staticmethod
    def _copy_file(src_path, dest_path):
        """Generic function to handle copying files."""
        if not os.path.exists(dest_path):
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)
            shutil.copy(src_path, dest_path)
        return dest_path

    def _convert_dates(self, array, year):
        """Convert every day in an array to days since Jan 1 1970."""
        # Loop through each position with data and convert
        ys, xs = np.where(array > 0)
        for y, x in zip(ys, xs):
            array[y, x] = self._convert_ordinal_to_unix_day(year, array[y, x])
        return array

    @staticmethod
    def _convert_ordinal_to_unix_day(year, ordinal_day):
        base = dt.datetime(1970, 1, 1)
        date = dt.datetime(year, 1, 1) + dt.timedelta(int(ordinal_day - 1))
        return (date - base).days

    @staticmethod
    def _convert_unix_day_to_calendar_date(unix_day):
        base = dt.datetime(1970, 1, 1)
        date = base + dt.timedelta(days=int(unix_day) - 1)
        return date.strftime("%Y-%m-%d")

    def _generate_land_cover_mosaic_dir(self, tile, year):
        return self.land_cover_dir.joinpath(tile,  str(year), "mosaics")

    def _generate_local_burn_hdf_dir(self, tile):
        return self.hdf_dir.joinpath(tile)

    def _generate_local_nc_path(self, tile):
        return self.nc_dir.joinpath(f"{tile}.nc")

    def _get_shape_files(self):
        """Get basic shapefiles needed for calculating statistics."""
        files_to_copy = {
            self._modis_sinusoidal_grid_shape_path: self.MODIS_SINUSOIDAL_PATH,
            self.conus_shape_path: self.CONUS_SHAPEFILE_PATH
        }
        for dest_path, source_path in files_to_copy.items():
            if not os.path.exists(dest_path):
                shutil.copy(source_path, dest_path)

    def _initialize_save_dirs(self):
        """Make all required project directories."""
        sdirs = [v for k, v in self.__dict__.items() if k.endswith("_dir")]
        for sdir in sdirs:
            os.makedirs(sdir, exist_ok=True)

    @staticmethod
    def _mode(vals) -> float:
        return max(set(list(vals)), key=list(vals).count)

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

    def _to_kms(self, p):
        return (p * self._res ** 2) / 1_000_000


class LPDAAC(Base):
    """Land Processes Distributed Active Archive Center access methods."""

    def __init__(self, project_directory, n_cores=0):
        """Initiate an LPDAAC object.

        Parameters
        ----------
        project_directory : str
            Path to firedpy project directory.
        """
        super().__init__(project_directory, n_cores)
        self._lp_daac_url = None
        self._date_regex = r"(?P<year>\d{4})\.(?P<month>\d{2})\.\
            (?P<day>\d{2})\/"

    def __repr__(self):
        """Return representation string for LPDAAC object."""
        name = self.__class__.__name__
        attrs = {}
        for key, attr in self.__dict__.items():
            if "data_frame" in key:
                attr = f"{type(attr)} {attr.shape}"  # Too big for preview
            if not key.startswith("_"):  # Avoid secrets/private attributes
                attrs[key] = attr
        address = hex(id(self))
        msgs = [f"\n   {k}='{v}'" for k, v in attrs.items()]
        msg = " ".join(msgs)
        return f"<{name} object at {address}> {msg}"

    def _download_files(self, download_requests):
        try:
            with Pool(self.n_cores - 1) as pool:
                mapper = pool.imap_unordered(
                    self._download_task,
                    download_requests
                )
                for _ in tqdm(mapper, total=len(download_requests)):
                    pass
        except Exception as pe:
            logger.error(f"Download failed: {pe}")
            try:
                for q in tqdm(download_requests, position=0, file=sys.stdout):
                    self._download_task(q)
            except Exception as e:
                template = "Download failed: error type {0}:\n{1!r}"
                message = template.format(type(e).__name__, e.args)
                logger.error(message)

    def _download_task(self, request: Tuple[str, str]):
        """Try request with EarthAccess, fallback to original method."""
        link = request[0]
        dest = request[1]
        if os.path.exists(dest):
            return

        # Try EarthAccess first if available
        if self._earthaccess.authenticated:
            try:
                success = self.download_file(link, dest)
                if success:
                    return
                else:
                    logger.error(
                        "EarthAccess download failed for "
                        f"{os.path.basename(dest)}, trying legacy method"
                    )
            except Exception as e:
                logger.error(
                    f"EarthAccess error: {e}, falling back to legacy method"
                )

        # Fallback to original urllib method
        try:
            pm = urllib.request.HTTPPasswordMgrWithDefaultRealm()
            pm.add_password(None, "https://urs.earthdata.nasa.gov",
                            self._username, self._password)

            cookie_jar = CookieJar()
            opener = urllib.request.build_opener(
                urllib.request.HTTPBasicAuthHandler(pm),
                urllib.request.HTTPCookieProcessor(cookie_jar)
            )
            urllib.request.install_opener(opener)

            myrequest = urllib.request.Request(link)
            response = urllib.request.urlopen(myrequest)

            with open(dest, "wb") as fd:
                while True:
                    chunk = response.read()
                    if chunk:
                        fd.write(chunk)
                    else:
                        break

        except Exception as e:
            logger.error(f"Download failed for {os.path.basename(dest)}: {e}")

    def _generate_local_hdf_path(self, year, remote_name):
        pass

    def _generate_tile(self, regex_group_dict: Dict[str, str]):
        horizontal_tile = regex_group_dict["horizontal_tile"]
        vertical_tile = regex_group_dict["vertical_tile"]
        return f"h{horizontal_tile}v{vertical_tile}"

    @staticmethod
    def get_all_available_tiles():
        with paramiko.SSHClient() as ssh_client:
            ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh_client.connect(hostname="fuoco.geog.umd.edu",
                               username="fire", password="burnt")
            logger.info("Connected to 'fuoco.geog.umd.edu' ...")
            # Open the connection to the SFTP
            with ssh_client.open_sftp() as sftp_client:
                sftp_client.chdir("/data/MODIS/C61/MCD64A1/HDF")
                return sftp_client.listdir()

    def _get_available_files(self, year_path: str, tiles: List[str] = None):
        url = urllib.parse.urljoin(self._lp_daac_url, year_path)
        request = requests.get(url)
        soup = BeautifulSoup(request.text, "html.parser")
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

    def _get_available_year_paths(self, start_year=None, end_year=None):
        # Get available years
        request = requests.get(self._lp_daac_url)
        soup = BeautifulSoup(request.text, "html.parser")
        year_paths = []
        for link in [link["href"] for link in soup.find_all("a", href=True)]:
            match = re.match(self._date_regex, link)
            if match is not None:
                file_year = int(match.groupdict().get("year"))
                if (start_year is None or file_year >= start_year) and (
                        end_year is None or file_year <= end_year):
                    year_paths.append(link)
        return year_paths

    @staticmethod
    def _verify_hdf_file(file_path):
        # Open the HDF file
        hdf_ds = gdal.Open(file_path)

        if hdf_ds is None:
            logger.warning(f"Failed to open {file_path}.")
            return False

        # List available sub-datasets (specific to HDF)
        sub_datasets = hdf_ds.GetSubDatasets()

        if not sub_datasets:
            logger.warning(f"No sub-datasets found in {file_path}.")
            return False

        return True


class BurnData(LPDAAC):
    """Methods for handling MODIS Burn data."""

    def __init__(self, project_directory, n_cores=0):
        """Initialize BurnData object.

        Parameters
        ----------
        project_directory : str | pathlib.PosixPath
            Target file output directory
        n_cores : int
            Number of CPU cores to use in multiprocessing. A value of 0 or
            None will use all available cores. Defaults to 0.
        """
        # Set these here, the LPDAAC will prompt if these aren't set yet
        super().__init__(project_directory=project_directory, n_cores=n_cores)

        self._lp_daac_url = "https://e4ftl01.cr.usgs.gov/MOTA/MCD64A1.061/"
        self._base_sftp_folder = os.path.join(
            "data", "MODIS", "C61", "MCD64A1", "HDF"
        )
        self._modis_template_path = os.path.join(
            project_directory, "rasters", "mosaic_template.tif"
        )
        self._record_start_year = 2000

    def __repr__(self):
        """Return representation string for BurnData object."""
        name = self.__class__.__name__
        attrs = {}
        for key, attr in self.__dict__.items():
            if "data_frame" in key:
                attr = f"{type(attr)} {attr.shape}"  # Too big for preview
            if not key.startswith("_"):  # Avoid secrets/private attributes
                attrs[key] = attr
        address = hex(id(self))
        msgs = []
        for key, value in attrs.items():
            if isinstance(value, (str, PosixPath)):
                msgs.append(f"\n   {key}='{value}'")
            else:
                msgs.append(f"\n   {key}={value}")
        msg = " ".join(msgs)
        return f"<{name} object at {address}> {msg}"

    def download_burn_data(self, tile, granules, download_dir):
        """Download a single tile's granules to a download diretory."""
        os.makedirs(download_dir, exist_ok=True)
        logger.info(f"Downloading to {tile} data to {download_dir}")
        try:
            downloaded_files = earthaccess.download(
                granules,
                download_dir
            )
            logger.info(
                f"Downloaded {tile} data ({len(downloaded_files)} files)"
            )

        except Exception as e:
            logger.error(f"Download failed for tile {tile}: {e}")

    def _create_requests(self, available_year_paths, tiles):
        """Create a list of requests for each target year.

        Parameters
        ----------
        available_year_paths : list
            List of paths associated with something....
        tiles : list
            List of tile strings (e.g., ['h12v04', 'h12v05'])

        Returns
        -------
        list : List of URL request strings.
        """
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
                        logger.info(f"Removing {local_file_path}")
                        os.remove(local_file_path)
                    else:
                        continue
                download_request = (
                    urllib.parse.urljoin(
                        urllib.parse.urljoin(self._lp_daac_url, year_path),
                        file
                    ),
                    local_file_path
                )
                download_requests.append(download_request)

        return download_requests

    def _extract_date_parts(self, path):
        fname = path.name
        match = re.match(self._file_regex, fname)
        if match:
            match_year = int(match.groupdict()["year"])
            match_day = int(match.groupdict()["ordinal_day"])
            return match_year, match_day
        return None

    def _generate_local_hdf_dir(self, tile):
        """Return local HDF directory path for given tile.

        Paramaters
        ----------
        tile : str
            Target MODIS tile (e.g., 'h12v04').

        Returns
        -------
        str : Path of HDF5 file for target MODIS tile.
        """
        return str(self.hdf_dir.joinpath(tile))

    def _generate_local_hdf_path(self, tile, hdf_name):
        """Return local HDF file path for given tile and HDF5 name.

        Paramaters
        ----------
        tile : str
            Target MODIS tile (e.g., 'h12v04').
        hdf_name : str
           Target MODIS HDF file name (e.g.,
           'MCD64A1.A2020245.h12v04.061.2021309112053.hdf').

        Returns
        -------
        str : Path of HDF5 file for target MODIS tile.
        """
        return str(self.hdf_dir.joinpath(tile, hdf_name))

    def get_burns(self, tiles=None, country=None, shape_file=None,
                  start_year=2000, end_year=2025):
        """Download MODIS burn data using EarthAccess.

        NOTE: One of tiles, country, or shapefile must be supplied. If mulitple
            of these parameters are provided, shape_file supercedes country and
            country supercedes tiles.

        Parameters
        ----------
        tiles : list
            List of MODIS tiles (e.g., ['h08v04', 'h09v04']).  Defaults to
            None.
        country : str
            The name of a country to use as a study area. Defaults to None.
        shape_file : str
            Path to a shapefile to use for the fire study area. Defaults to
            None.
        start_year : int
            The first year of fire events. Defaults to 2000.
        end_year : int
            The last year of fire events. Defaults to 2025.
        """
        start_date = f"{start_year}-01-01"
        end_date = f"{end_year}-12-31"
        logger.info(f"Getting burn data using EarthAccess for tiles: {tiles}")
        logger.info(f"Date range: {start_date} to {end_date}")

        try:
            # Search for granules
            logger.info(f"Searching for MCD64A1 granules: {start_date} - "
                        f"{end_date}")
            granule_dict = self._get_granules(tiles, country, shape_file,
                                              start_date, end_date)
            tiles = list(granule_dict)

            # Collect each tile granules and download directory
            tile_granules = {}
            download_dirs = {}
            for tile, granules in granule_dict.items():
                logger.info(f"Processing tile: {tile}")

                # Check if any granules were found
                n_granules = len(granules)
                logger.info(f"Found {n_granules} granules for tile {tile}")
                if not granules:
                    logger.info(f"No data found for tile {tile}")
                    continue

                # Collect results
                download_dir = self.hdf_dir.joinpath(tile)
                download_dirs[tile] = download_dir
                tile_granules[tile] = granules

            # Download granules (3 Core limit with Earthdata)
            n_cores = min(3, self.n_cores)
            if n_cores > 1:
                with ThreadPoolExecutor(n_cores) as pool:
                    jobs = []
                    for tile in tiles:
                        download_dir = download_dirs[tile]
                        granules = tile_granules[tile]
                        job = pool.submit(
                            self.download_burn_data,
                            tile,
                            granules,
                            download_dir
                        )
                        jobs.append(job)
                    for job in as_completed(jobs):
                        job.result()
            else:
                for tile in tiles:
                    download_dir = download_dirs[tile]
                    granules = tile_granules[tile]
                    self.download_burn_data(
                        tile,
                        granules,
                        download_dir
                    )

            logger.info("EarthAccess burn data download completed!")

        except Exception as e:
            logger.error(f"EarthAccess download failed: {e}")
            traceback.print_exc()
            raise RuntimeError(
                f"Failed to get burn data with EarthAccess: {e}"
            )

        # Convert to NetCDF
        self._write_ncs(tiles)
        logger.info(f"Created NetCDF for tile(s) {tiles}")

        return tiles

    def _get_granules(self, tiles, country, shape_file, start_date, end_date):
        """Get Earth Access data granules using tile names."""
        # Get tile list
        if shape_file:
            tiles = shape_to_tiles(shape_file)
        elif country:
            tiles = country_to_tiles(country)
        else:
            if not tiles:
                raise KeyError(
                    "One of `tiles`, `country`, or `shape_file` parameters"
                    "must be supplied."
                )

        # Let's see if using points is quicker
        points = tiles_to_points(tiles)

        # We can use up to three cores as per Earthdata rules
        n_cores = min(3, self.n_cores)
        granules = {}
        if n_cores > 1:
            with ThreadPoolExecutor(n_cores) as pool:
                jobs = {}
                for tile, point in points.items():
                    job = pool.submit(
                        earthaccess.search_data,
                        short_name="MCD64A1",
                        version="061",
                        temporal=(start_date, end_date),
                        point=point
                    )
                    jobs[tile] = job
                for tile, job in tqdm(jobs.items(), total=len(jobs)):
                    granule = job.result()
                    granules[tile] = granule
        else:
            for tile, point in tqdm(points.items(), total=len(points)):
                granule = earthaccess.search_data(
                    short_name="MCD64A1",
                    version="061",
                    temporal=(start_date, end_date),
                    point=point
                )
                granules[tile] = granule

        return granules

    def _get_search_kwargs(self, tile, start_date, end_date):
        kwargs = dict(
            short_name="MCD64A1",
            version="061",
            temporal=(start_date, end_date),
            granule_name=f"*{tile}*"
        )
        return kwargs

    def _write_modis_template_file(self):
        # Merge one year into a reference mosaic
        logger.info(f"Creating reference mosaic raster...")
        folders = glob(os.path.join(self.hdf_dir, "*"))
        file_groups = [glob(os.path.join(f, "*hdf")) for f in folders]
        for f in file_groups:
            f.sort()
        files = [f[0] for f in file_groups]
        dss = [rasterio.open(f).subdatasets[0] for f in files]
        tiles = [rasterio.open(d) for d in dss]
        mosaic, transform = merge(tiles)
        crs = tiles[0].meta.copy()
        crs.update(
            {
                "driver": "GTIFF",
                "height": mosaic.shape[1],
                "width": mosaic.shape[2],
                "transform": transform
            }
        )

        dst = self._modis_template_path
        logger.info(f"Writing reference mosaic raster to {dst}.")
        with rasterio.open(dst, "w+", **crs) as r:
            r.write(mosaic)

    @staticmethod
    def _verify_hdf_file(file_path):
        # Open the HDF file
        hdf_ds = gdal.Open(file_path)

        if hdf_ds is None:
            logger.warning(f"Failed to open {file_path} with GDAL.")
            return False

        # List available sub-datasets (specific to HDF)
        sub_datasets = hdf_ds.GetSubDatasets()

        if not sub_datasets:
            logger.warning(f"No sub-datasets found in {file_path}.")
            return False

        return True

    def _write_ncs(self, tiles):
        """Convert a list of tiles associated with MODIS HDF4 to a NetCDF file.

        Parameters
        ----------
        tiles : list[str]
            A list of strings representing MODIS tile IDs.
        """
        # Build the netcdfs here
        fill_value = -9999  # Default fill value in original file is -1
        logger.info(f"Building netcdf files for tiles {tiles}")
        for tile in tqdm(tiles):

            # Build/find file paths needed for this operation
            nc_file_path = self._generate_local_nc_path(tile)
            hdf_dir = self.hdf_dir.joinpath(tile)
            if nc_file_path.exists():
                continue
            paths = []
            for path in list(hdf_dir.glob("*.hdf")):
                if self._extract_date_parts(path):
                    paths.append(path)
            files = sorted(paths, key=self._extract_date_parts)

            if not files:
                msg = f"No HDF4 files for tile {tile} in {hdf_dir}"
                logger.error(msg)
                raise OSError(msg)

            try:
                # Use a sample to get geography information and geometries
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
                nco = Dataset(
                    nc_file_path,
                    mode="w",
                    format="NETCDF4",
                    clobber=True
                )

                # Dimensions
                nco.createDimension("y", ny)
                nco.createDimension("x", nx)
                nco.createDimension("time", None)

                # Variables
                y = nco.createVariable("y", np.float64, ("y",))
                x = nco.createVariable("x", np.float64, ("x",))
                times = nco.createVariable("time", np.int16, ("time",))
                variable = nco.createVariable(
                    "value",
                    np.int16,
                    ("time", "y", "x"),
                    fill_value=fill_value,
                    zlib=True
                )
                variable.standard_name = "day"
                variable.long_name = "Burn Days"

                # Appending the CRS information
                # Check "https://cf-trac.llnl.gov/trac/ticket/77"
                crs = nco.createVariable("crs", "c")
                variable.setncattr("grid_mapping", "crs")
                crs.spatial_ref = proj
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
                days = []
                for file in files:
                    p = self._extract_date_parts(file)
                    day = self._convert_ordinal_to_unix_day(p[0], p[1])
                    days.append(day)
                days = np.array(days)

                # Write dimension data
                x[:] = xs
                y[:] = ys
                times[:] = days

                # One file a time, write the arrays
                for tile_index, f in enumerate(files):
                    match = re.match(self._file_regex, os.path.basename(f))
                    if match is None:
                        continue

                    regex_group_dict = match.groupdict()

                    try:
                        ds = gdal.Open(f).GetSubDatasets()[0][0]
                        hdf = gdal.Open(ds)
                    except Exception as e:
                        logger.error(
                            f"Could not open {f} for building ncdf: {e}"
                        )
                        variable[tile_index, :, :] = np.full(
                            (ny, nx),
                            fill_value
                        )
                        continue

                    data = hdf.GetRasterBand(1)
                    array = data.ReadAsArray()
                    nulls = np.where(array < 0)

                    year = int(regex_group_dict["year"])
                    array = self._convert_dates(array, year)
                    array[nulls] = fill_value

                    if array.shape == (ny, nx):
                        variable[tile_index, :, :] = array
                    else:
                        logger.warning(
                            f"{f}: failed, had wrong dimensions, "
                            "inserting a blank array in its place."
                        )
                        variable[tile_index, :, :] = np.full(
                            (ny, nx),
                            fill_value
                        )

                # Done
                nco.close()

            except Exception as e:
                # Log the error and move on to the next tile
                logger.error(f"Error processing tile {tile}: {str(e)}")
                raise OSError(f"Error processing tile {tile}: {str(e)}")


class LandCover(Base):
    """EarthAccess-based Land Cover data access for firedpy.

    TODO: Double check that the land cover and fired datasets line up
        correctly. There may be discrepancies in the years.
    """

    def __init__(
            self,
            project_directory,
            n_cores=0
    ):
        super().__init__(project_directory, n_cores)

        # Smart mapping in cases when land cover tiles differ from burn tiles
        self._tile_mapping = {
            "h09v04": ["h09v03", "h08v04"],  # Try nearby tiles
            "h09v05": ["h09v06", "h08v04"],
            "h10v04": ["h10v02", "h09v03"],
            "h10v05": ["h10v06", "h09v06"]
        }

    def __repr__(self):
        """Return representation string for a LandCover object."""
        name = self.__class__.__name__
        attrs = {}
        for key, attr in self.__dict__.items():
            if "data_frame" in key:
                attr = f"{type(attr)} {attr.shape}"  # Too big for preview
            if not key.startswith("_"):  # Avoid secrets/private attributes
                attrs[key] = attr
        address = hex(id(self))
        msgs = [f"\n   {k}='{v}'" for k, v in attrs.items()]
        msg = " ".join(msgs)
        return f"<{name} object at {address}> {msg}"

    def add_land_cover_attributes(self, gdf, tiles, land_cover_type=0):
        """Add landcover attributes to a firedpy GeoDataFrame.

        Parameters
        ----------
        gdf : geopandas.geodataframe.GeoDataFrame
            A geodata frame of fire events.
        tiles : list[str]
            A list of MODIS of grid tile IDs.
        land_cover_type : int | LandCoverType
            Include land cover as an attribute, provide a number corresponding
            with a MODIS/Terra+Aqua Land Cover (MCD12Q1) category followed with
            username:password of your NASA's Earthdata service account.
            Available land cover categories:

                1: IGBP global vegetation classification scheme
                2: University of Maryland (UMD) scheme
                3: MODIS-derived LAI/fPAR scheme

            If you do not have an account register at
            https://urs.earthdata.nasa.gov/home. Defaults to 0 or
            LandCoverType.NONE.

        Returns
        -------
        geopandas.geodataframe.GeoDataFrame : The geodataframe with landcover
            attributes added.
        """
        if gdf.shape[0] == 0:
            logger.warning(
                f"No fire events found in tiles {tiles}, not adding land cover"
                " attributes"
            )
            return gdf

        # This will use a config file or a prompt if one isn't available
        auth = earthaccess.login(strategy="all")
        if not auth:
            msg = "EarthAccess authentication failed"
            logger.error(msg)
            raise RuntimeError(msg)
        else:
            logger.info("EarthAccess authentication successful.")

        # Match possible land cover type argument types
        logger.info("Adding land cover attributes...")
        if isinstance(land_cover_type, str):
            land_cover_type = int(land_cover_type)
        if isinstance(land_cover_type, LandCoverType):
            land_cover_type = land_cover_type.value

        # We'll need to specify which type of land_cover
        lc_descriptions = {
            1: "IGBP global vegetation classification scheme",
            2: "University of Maryland (UMD) scheme",
            3: "MODIS-derived LAI/fPAR scheme",
            4: "MODIS-derived Net Primary Production (NPP) scheme",
            5: "Plant Functional Type (PFT) scheme."
        }

        # Rasterio point querier (why define here?)
        def point_query(row):
            x = row["x"]
            y = row["y"]
            try:
                val = [val for val in lc.sample([(x, y)])][0][0]
            except Exception:
                val = np.nan
            return val

        # Get the range of burn years
        burn_years = list(gdf["ig_year"].unique())
        burn_years.sort()

        # This works faster when split by year and the pointer is outside
        # This is also not the best way
        sgdfs = []
        for tile in tiles:
            for year in tqdm(burn_years, position=0, file=sys.stdout):
                # Get the land cover geotiff directory for this tile
                mosaic_dir = self.land_cover_dir.joinpath(
                    tile, str(year), "mosaics"
                )
                if not mosaic_dir.exists():
                    logger.warning(
                        f"No land cover data for {tile} in year {year}"
                    )
                    continue

                # I don't understand why were taking this step below  # <------ Address this
                # lc_files = []
                # lc_years = []
                # for fname in os.listdir(mosaic_dir):
                #     # How necessary is the re match?
                #     # if re.match(self._lc_mosaic_re, fname):
                #     fpath = os.path.join(mosaic_dir, fname)
                #     lc_files.append(fpath)
                # lc_files = sorted(lc_files)

                # # Group by year?
                # for fpath in lc_files:
                #     fname = os.path.basename(fpath)
                #     # Don't we know what year it is?
                #     # year_match = re.match(self._lc_mosaic_re, fname)
                #     # year_dict = year_match.groupdict()
                #     year = int(year_dict["year"])
                #     lc_years.append(year)
                # lc_files = {lc_years[i]: f for i, f in enumerate(lc_files)}

                # Now set year one back for land_cover
                # year = year - 1

                # Use previous year's land cover
                # if year < min(lc_years):
                #     year = min(lc_years)
                # elif year > max(lc_years):
                #     year = max(lc_years)

                # lc_file = lc_files[year]

                # Collect all files in the mosaic directory (all? there's 1)
                lc_file = list(mosaic_dir.glob("*tif"))
                if not lc_file:
                    continue
                else:
                    lc_file = lc_file[0]

                # Create a sub data frame for this year
                sgdf = gdf[gdf["ig_year"] == year].copy()

                # Open the land cover raster and add attributes to sub df
                lc = rasterio.open(lc_file)
                sgdf.loc[:, "lc_code"] = sgdf.apply(point_query, axis=1)
                sgdf = sgdf[sgdf["lc_code"] != 255]  # Out-of-tile points
                idgrp = sgdf.groupby("id")
                sgdf.loc[:, "lc_mode"] = idgrp["lc_code"].transform(self._mode)
                sgdfs.append(sgdf)

        gdf = pd.concat(sgdfs)
        gdf = gdf.reset_index(drop=True)

        # Add in the class description from land_cover tables
        land_cover_path = self._copy_land_cover_ref(land_cover_type)
        lc_table = pd.read_csv(land_cover_path)
        gdf = pd.merge(
            left=gdf,
            right=lc_table,
            how="left",
            left_on="lc_mode",
            right_on="Value"
        )
        gdf = gdf.drop("Value", axis=1)
        gdf.loc[:, "lc_type"] = lc_descriptions[land_cover_type]
        gdf.rename({"lc_description": "lc_desc"}, inplace=True, axis="columns")

        return gdf

    def _copy_land_cover_ref(self, land_cover_type: int) -> str:
        lookup = DATA_DIR.joinpath(
            "land_cover",
            f"MCD12Q1_LegendDesc_Type{land_cover_type}.csv"
        )
        land_cover_out_dir_path = self.project_directory.joinpath(
            "tables", "land_cover",
            f"MCD12Q1_LegendDesc_Type{land_cover_type}.csv"
        )
        return self._copy_file(lookup, land_cover_out_dir_path)

    def _copy_wwf_file(self) -> str:
        lookup = DATA_DIR.joinpath('world_eco_regions', 'wwf_terr_ecos.gpkg')
        wwf_out_dir_path = os.path.join(
            self.project_directory, 'shapefiles', 'eco_region',
            'wwf_terr_ecos.gpkg'
        )
        return self._copy_file(lookup, wwf_out_dir_path)

    def _find_available_tiles_for_region(self, requested_tiles):
        """Find available land cover tiles that cover the requested region."""
        try:
            # Get all available tiles for a sample year
            granules = earthaccess.search_data(
                short_name="MCD12Q1",
                version="061",
                temporal=("2020-01-01", "2020-12-31"),  # Why 2020?
                count=500
            )

            available_tiles = set()
            for granule in granules:
                name = granule.get("meta", {}).get("native-id", "")
                tile_match = re.search(r"\.h(\d{2})v(\d{2})\.", name)
                if tile_match:
                    h, v = tile_match.groups()
                    tile = f"h{h}v{v}"
                    available_tiles.add(tile)

            # For each requested tile, find available alternatives
            tiles_to_use = set()
            for tile in requested_tiles:
                if tile in available_tiles:
                    tiles_to_use.add(tile)
                    logger.info(f"   {tile} - Available directly")
                elif tile in self._tile_mapping:
                    for alt_tile in self._tile_mapping[tile]:
                        if alt_tile in available_tiles:
                            tiles_to_use.add(alt_tile)
                            logger.info(
                                f"   {tile} -> using {alt_tile} (nearby "
                                "coverage)"
                            )
                            break
                    else:
                        logger.warning(
                            f"   {tile} - No suitable alternative found"
                        )
                else:
                    logger.warning(f"   {tile} - No mapping defined")

            return list(tiles_to_use)

        except Exception as e:
            logger.error(f"Error finding available tiles: {e}")
            return []

    def _generate_land_cover_mosaic_dir(self, tile: str, year: str) -> str:
        return os.path.join(self.land_cover_dir, tile, str(year), "mosaics")

    def _generate_local_hdf_path(self, tile, year, remote_name):
        return os.path.join(self.land_cover_dir, tile, year, remote_name)

    def _generate_local_hdf_dir(self, tile: str, year: str) -> str:
        return os.path.join(self.land_cover_dir, tile, year)

    def get_land_cover(self, gdf, tiles, land_cover_type=1):
        """Download and process land cover data with EarthAccess.

        Parameters
        ----------
        gdf : geopandas.geodataframe.GeoDataFrame
            A geodata frame of fire events.
        tiles : list
            List of MODIS tiles (e.g., ['h08v04', 'h09v04']).
        land_cover_type : int
            Include land cover as an attribute, provide a number corresponding
            with a MODIS/Terra+Aqua Land Cover (MCD12Q1) category followed with
            username:password of your NASA's Earthdata service account.
            Available land cover categories:

                1: IGBP global vegetation classification scheme
                2: University of Maryland (UMD) scheme
                3: MODIS-derived LAI/fPAR scheme
                4: Annual BIOME-Biogeochemical Cycles (BGC)
                5: Annual Plant Functional Types (PFT)

            If you do not have an account register at
            https://urs.earthdata.nasa.gov/home. Defaults to 1.
        """
        if not self._earthaccess.authenticated:
            logger.warning("EarthAccess not authenticated for land cover.")
            return

        if tiles is None:
            logger.warning("No tiles specified for land cover.")
            return

        # Find available tiles that can cover our region
        logger.info(
            f"Getting land cover data using EarthAccess for region: {tiles}"
        )
        logger.info("Finding available land cover tiles for region...")
        available_tiles = self._find_available_tiles_for_region(tiles)

        if not available_tiles:
            logger.warning("No land cover tiles available for this region.")
            return

        try:
            # Get available years (this is fixed?)
            available_years = gdf["ig_year"].unique()
            for tile in available_tiles:
                logger.info(f"\nProcessing land cover for tile: {tile}")

                for year in available_years:
                    logger.info(f"   Processing year: {year}")
                    year = str(year)

                    # Check if mosaic already exists
                    lc_dir = self.land_cover_dir
                    mosaic_dir = lc_dir.joinpath(tile, year, "mosaics")
                    tag = f"{tile}_{land_cover_type}_{year}"
                    output_file = f"lc_mosaic_{tag}.tif"
                    mosaic_path = os.path.join(mosaic_dir, output_file)
                    os.makedirs(mosaic_dir, exist_ok=True)
                    if os.path.exists(mosaic_path):
                        logger.info(f"   Mosaic already exists: {output_file}")
                        continue

                    # Search for granules for this year and tile
                    logger.info(f"Searching for {tile} data in {year}...")
                    granules = earthaccess.search_data(
                        short_name="MCD12Q1",
                        version="061",
                        temporal=(f"{year}-01-01", f"{year}-12-31")
                    )

                    # Filter for this specific tile
                    tile_granules = []
                    for granule in granules:
                        id = "native-id"
                        granule_name = granule.get("meta", {}).get(id, "")
                        if tile in granule_name:
                            tile_granules.append(granule)
                    if not tile_granules:
                        msg = f"No land cover data for {tile} in year {year}."
                        logger.warning(msg)
                        continue
                    msg = f"Found {len(tile_granules)} granules, downloading.."
                    logger.info(msg)

                    # Download granules
                    download_dir = self._generate_local_hdf_dir(tile, year)
                    os.makedirs(download_dir, exist_ok=True)
                    try:
                        downloaded_files = earthaccess.download(
                            granules=tile_granules,
                            local_path=download_dir
                        )
                        msg = f"Downloaded {len(downloaded_files)} files."
                        logger.info(msg)

                        # Create mosaic (simplified version)
                        if downloaded_files:
                            msg = f"Land cover data ready for {tile} {year}."
                            logger.info(msg)

                            # This only creates a single netcdf per year and
                            # tile. The way the filesytem and loop here is
                            # organized, that's all we can do here. Later,
                            # the method that adds this attribute to the
                            # dataframe assumes it is a full study area mosaic
                            # so we are simply dropping missed event points,
                            # which isn't ideal. Though it works, it doesn't
                            # seem like this was the original intention.
                            self.mosaic_landuse(
                                downloaded_files,
                                mosaic_path,
                                land_cover_type
                            )

                    except Exception as e:
                        msg = f"Processing failed for {tile} {year}: {e}"
                        logger.error(msg)
                        continue

            logger.info("EarthAccess land cover processing completed.")

        except Exception as e:
            logger.error(f"EarthAccess land cover failed: {e}")
            logger.info("Continuing without land cover data...")

    def mosaic_landuse(self, downloaded_files, mosaic_path, land_cover_type):
        """Merge a list of MODIS HDF files into one GeoTiff.

        TODO: Will there ever be multiple downloaded files from
            LandCover.get_land_cover? I was assuming there would be since
            we are labeling the final output tiff "mosaic" but it looks like
            the earthaccess method that downloads the originals takes in
            one tile at a time.

        Parameters
        ----------
        downloaded_files : list[str]
            A list of paths to MODIS Land Cover HDF4 files.
        mosaic_path : str | pathlib.PosixPath
            Target file path for the GeoTiff.
        land_cover_type : int | firedpy.enums.LandCoverType
            Include land cover as an attribute, provide a number corresponding
            with a MODIS/Terra+Aqua Land Cover (MCD12Q1) category followed with
            username:password of your NASA's Earthdata service account.
            Available land cover categories:

                1: IGBP global vegetation classification scheme
                2: University of Maryland (UMD) scheme
                3: MODIS-derived LAI/fPAR scheme

            If you do not have an account register at
            https://urs.earthdata.nasa.gov/home. Defaults to 1.
        """
        # Use the target file's parent directory for the temporary files
        mosaic_dir = Path(os.path.dirname(mosaic_path))
        tmp_files = []
        for i, file in enumerate(downloaded_files):
            # Get a variable name pattern
            pattern = f"LC_Type{land_cover_type}"

            # Create temporary file
            dst = mosaic_dir.joinpath(f"tile_{i}.tif")
            hdf4_to_geotiff(file, dst, pattern=pattern)
            tmp_files.append(dst)

        # Mosaic these together and save to file
        if len(tmp_files) > 1:
            shutil.rmtree(tmp_files)
            raise NotImplementedError(
                "Merging multiple land use layers not implemented yet: "
                f"{downloaded_files}"
            )

        else:
            shutil.move(tmp_files[0], mosaic_path)


class EcoRegion(Base):
    """Methods for managing Ecoregion data."""

    def __init__(self, project_directory):
        """Iniitialize an EcoRegion object.

        Parameters
        ----------
        project_directory : str | pathlib.PosixPath
            Path to firedpy output directory.
        """
        super().__init__(project_directory)
        self._eco_region_ftp_url = (
            "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/"
            "Ecoregions/cec_na/NA_CEC_Eco_Level3.zip"
        )
        self._eco_region_raster_path = os.path.join(
            self.eco_region_raster_dir,
            "NA_CEC_Eco_Level3_modis.tif"
        )
        self._ref_cols = [
            "NA_L3CODE", "NA_L3NAME", "NA_L2CODE", "NA_L2NAME", "NA_L1CODE",
            "NA_L1NAME", "NA_L3KEY", "NA_L2KEY", "NA_L1KEY"
        ]
        self.eco_region_data_frame = None

    def __repr__(self):
        """Return representation string for an EcoRegion object."""
        name = self.__class__.__name__
        attrs = {}
        for key, attr in self.__dict__.items():
            if "data_frame" in key:
                if attr is not None:
                    attr = f"{type(attr)} {attr.shape}"  # Too big for preview
            if not key.startswith("_"):  # Avoid secrets/private attributes
                attrs[key] = attr
        address = hex(id(self))
        msgs = [f"\n   {k}='{v}'" for k, v in attrs.items()]
        msg = " ".join(msgs)
        return f"<{name} object at {address}> {msg}"

    def add_eco_region_attributes(
            self,
            gdf,
            eco_region_type="na",
            eco_region_level=None
    ):
        """Add eco region attributes to a fire event geodataframe.

        Parameters
        ----------
        gdf : geopandas.geodataframe.GeoDataFrame
        eco_region_type : str
            Specify the ecoregion type as either 'world' or 'na':

                'world' = World Terrestrial Ecoregions (World Wildlife Fund)
                'na' = North American ecoregions (Omernick, 1987)

            The most common (modal) ecoregion across the event is used.

            Further, to associate each event with North American ecoregions
            (Omernick, 1987) you may provide a number corresponding to an
            ecoregion level. Ecoregions are retrieved from www.epa.gov and
            levels I through IV are available. Levels I and II were developed
            by the North American Commission for Environmental Cooperation.
            Levels III and IV were developed by the United States Environmental
            Protection Agency. For events with more than one ecoregion, the
            most common value will be used. Defaults to None.
        eco_region_level : int
            The desired Ecoregions level from the North American Commission for
            Environmental Cooperation (CEC). Levels 1 to 3 are available, with
            level 1 representing the broadest scale and level III representing
            the most detailed. Defaults to 1.

        Returns
        -------
        geopandas.geodataframe.GeoDataFrame: A copy of the geodataframe with
            eco region attributes added.
        """
        if gdf.shape[0] == 0:
            msg = "No fire events found, not adding eco-region attributes."
            logger.error(msg)
            return gdf

        logger.info("Adding eco-region attributes ...")
        eco_region_type = EcoRegionType(eco_region_type)
        if eco_region_level or (eco_region_type == EcoRegionType.NA):
            gdf = self.add_attributes_from_na_cec(gdf, eco_region_level)
        else:
            gdf = self.add_attributes_from_wwf(gdf)
        return gdf

    def add_attributes_from_na_cec(self, gdf, eco_region_level):
        # Different levels have different sources
        institution = "(NA-Commission for Environmental Cooperation)"
        eco_types = {
            "NA_L1CODE": f"Level I Ecoregions {institution}",
            "NA_L2CODE": f"Level II Ecoregions {institution}",
            "NA_L3CODE": f"Level III Ecoregions {institution}"
        }

        # Read in the Level File (contains every level) and reference table
        shp_path = self._copy_cec_file()
        eco = gpd.read_file(shp_path)
        eco.to_crs(gdf.crs, inplace=True)

        # Filter for selected level (level III defaults to US-EPA version)
        if not eco_region_level:
            eco_region_level = 3
        eco_code = [c for c in eco if str(eco_region_level) in
                    c and "CODE" in c]

        if len(eco_code) > 1:
            eco_code = [c for c in eco_code if "NA" in c][0]
        else:
            eco_code = eco_code[0]

        logger.info(f"Selected ecoregion code: {eco_code}")

        # Find modal eco-region for each event id
        eco = eco[[eco_code, "geometry"]]
        gdf = gpd.sjoin(gdf, eco, how="left", predicate="within")
        gdf = gdf.reset_index(drop=True)
        gdf["eco_mode"] = gdf.groupby("id")[eco_code].transform(self._mode)

        # Add in the type of eco-region
        gdf["eco_type"] = eco_types[eco_code]

        # Add in the name of the modal ecoregion
        eco_ref = pd.read_csv(self.eco_region_csv_path)
        eco_name = eco_code.replace("CODE", "NAME")
        eco_df = eco_ref[[eco_code, eco_name]].drop_duplicates()
        eco_map = dict(zip(eco_df[eco_code], eco_df[eco_name]))
        gdf["eco_name"] = gdf[eco_code].map(eco_map)

        # Clean up column names
        gdf = gdf.drop("index_right", axis=1)
        gdf = gdf.drop(eco_code, axis=1)

        return gdf

    def add_attributes_from_wwf(self, gdf):
        # Read in the world ecoregions from WWF
        eco_path = self._copy_wwf_file()
        eco = gpd.read_file(eco_path)
        eco.to_crs(gdf.crs, inplace=True)

        # Find modal eco region for each event id
        eco = eco[["ECO_NUM", "ECO_NAME", "geometry"]]
        gdf = gpd.sjoin(gdf, eco, how="left", predicate="within")
        gdf = gdf.reset_index(drop=True)

        gdf["eco_mode"] = gdf.groupby("id")["ECO_NUM"].transform(self._mode)
        gdf["eco_name"] = gdf["ECO_NAME"]
        gdf["eco_type"] = "WWF Terrestrial Ecoregions of the World"

        gdf = gdf.drop("index_right", axis=1)
        gdf = gdf.drop("ECO_NAME", axis=1)
        gdf = gdf.drop("ECO_NUM", axis=1)

        return gdf

    def create_eco_region_raster(self, tiles: List[str]):
        """Create an EcoRegion raster.

        Parameters
        ----------
        tiles: list
            A list of strings representing the target MODIS tiles in which to
            create the Ecoregion raster, e.g., ["", ""]

        Returns
        -------
        """
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
            burn_dir = os.path.join(self.hdf_dir, tile)
            if not os.path.exists(burn_dir):
                raise FileNotFoundError(
                    f"Burn HDFs do not exist for tile {tile} at {burn_dir}. "
                    "Use the BurnData class to download this data before "
                    "rasterizing the eco region file"
                )

            # Find the matching file
            files = []
            for f in glob(os.path.join(burn_dir, "*")):
                if re.match(self._file_regex, os.path.basename(f)) is not None:
                    files.append(os.path.join(burn_dir, f))
            file = files[0]
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
        self._rasterize_vector_data(
            self.eco_region_data_frame,
            self._eco_region_raster_path,
            attribute,
            xres,
            projection,
            extent
        )

    def get_eco_region(self):
        """Download Ecoregion shapefile to the project directory.

        NOTE: This currently only downloads the EPA Omernick Ecoregions North
            for North American, though we tell the user they could have this
            or the World Terrestrial Ecoregions (World Wildlife Fund).
        """
        # Download the file from source
        eco = self._read_eco_region_file()

        # Create a reference table for ecoregions
        eco_ref = eco[self._ref_cols].drop_duplicates()
        eco_ref = eco_ref.map(self._normalize_string)
        eco_ref.to_csv(self.eco_region_csv_path, index=False)
        self.eco_region_data_frame = eco_ref

    @staticmethod
    def _normalize_string(string):
        def capitalize_special(s):
            """Handle capitalization for special characters within a string."""
            if "/" in s:
                segments = []
                for segment in s.split("/"):
                    if segment.upper() != "USA":
                        segment = segment.capitalize()
                    else:
                        segment = segment.upper()
                    segments.append(segment)
                s = "/".join(segments)
            if "-" in s:
                segments = []
                for segment in s.split("-"):
                    if segment.upper() != "USA":
                        segment = segment.title()
                    else:
                        segment = segment.upper()
                    segments.append(segment)
                s = "-".join(segments)
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

    def _read_eco_region_file(self):
        """Read Ecoregion file from package data or EPA FTP site.

        NOTE: Update 02/2021: EPA FTP site is glitchy, using local file in case
        download fails.

        Returns
        -------
        gpd.
        """
        # Check if the file exists, try to use a new EPA file
        if not os.path.exists(self.eco_region_shape_path):
            try:
                eco = gpd.read_file(self._eco_region_ftp_url)
                eco = eco.to_crs("epsg:5070")
                eco.to_file(self.eco_region_shape_path)
            except Exception as e:
                msg = f"Download from EPA FTP site, using local file {e}"
                logger.info(msg)
                shutil.copytree(
                    self.project_eco_region_dir,
                    self.eco_region_shape_dir
                )

        # Read in the file downloaded or copied to the output directory
        df = gpd.read_file(self.eco_region_shape_path)

        return df
