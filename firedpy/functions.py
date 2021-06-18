# -*- coding: utf-8 -*-
from bs4 import BeautifulSoup
from collections import OrderedDict
import datetime as dt
import gc
import geopandas as gpd
from getpass import getpass
from glob import glob
from io import BytesIO
from multiprocessing import cpu_count, Pool
from netCDF4 import Dataset
import numpy as np
import os
import pandas as pd
import pycurl
import rasterio
from rasterio import logging
from rasterio.merge import merge
import xarray as xr
from shapely.geometry import Point, Polygon, MultiPolygon
import sys
from tqdm import tqdm
import requests
import warnings
import paramiko

# The python gdal issue (matching system gdal version)
try:
    from osgeo import gdal, ogr, osr
except ImportError:
    raise ImportError(""" Unfortunately, you still need to install GDAL for
                      Python. Try pip install `pygdal==version` where the
                      version matches the first three digits of the output from
                      the command `gdalinfo --version`. To see available pygdal
                      versions run `pip install pygdal== '
                      """)

# MODIS CRS retrieved from a single HDF file
outCRS = '''PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",
DATUM["Not specified (based on custom spheroid)",
SPHEROID["Custom spheroid",6371007.181,0]],
PRIMEM["Greenwich",0],
UNIT["degree",0.0174532925199433,
AUTHORITY["EPSG","9122"]]],
PROJECTION["Sinusoidal"],
PARAMETER["longitude_of_center",0],
PARAMETER["false_easting",0],
PARAMETER["false_northing",0],
UNIT["Meter",1],
AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''

# Suppress rasterio errors for now
log = logging.getLogger()
log.addFilter(rasterio.errors.NotGeoreferencedWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
pd.options.mode.chained_assignment = None


# Functions
def convertDates(array, year):
    """Convert everyday in an array to days since Jan 1 1970"""

    def convertDate(julien_day, year):
        base = dt.datetime(1970, 1, 1)
        date = dt.datetime(year, 1, 1) + dt.timedelta(int(julien_day))
        days = date - base
        return days.days

    # Loop through each position with data and convert
    locs = np.where(array > 0)
    ys = locs[0]
    xs = locs[1]
    locs = [[ys[i], xs[i]] for i in range(len(xs))]
    for loc in locs:
        y = loc[0]
        x = loc[1]
        array[y, x] = convertDate(array[y, x], year)

    return array


def dateRange(perimeter):
    """Converts days in a perimeter object since Jan 1 1970 to date strings"""

    if len(perimeter.coords) > 0:
        base = dt.datetime(1970, 1, 1)
        days = [p[2] for p in perimeter.coords]
        day1 = (base + dt.timedelta(days=int(min(days)))).strftime("%Y-%m-%d")
    else:
        day1 = "N/A"
    return day1


def edgeCheck(yedges, xedges, coord, sp_buffer):
    """Identify edge cases to make merging events quicker later"""

    y = coord[0]
    x = coord[1]
    if y in yedges:
        edge = True
    elif x in xedges:
        edge = True
    else:
        edge = False
    return edge


def flttn(lst):
    """Just a quick way to flatten lists of lists"""
    lst = [l for sl in lst for l in sl]
    return lst


# def maxGrowthDate(x):
#     dates = x["date"].to_numpy()
#     pixels = x["pixels"].to_numpy()
#     loc = np.where(pixels == np.max(pixels))[0]
#     d = np.unique(dates[loc])[0]
#     # if len(d) > 1:
#     #     d = ", ".join(d)
#     # else:
#     #     d = d[0]
#     return d


def maxGrowthDate(x):
    dates = x["date"].to_numpy()
    pixels = x["pixels"].to_numpy()
    loc = np.where(pixels == np.max(pixels))[0]
    d = np.unique(dates[loc])[0]
    # if len(d) > 1:
    #     d = ", ".join(d)
    # else:
    #     d = d[0]
    return d


def mergeChecker(new_coords, full_list, temporal_param, radius):
    """
    This uses a radius for the spatial window as opposed to a square and is not
    currently being used to merge events.
    """
    t1 = np.min([c[2] for c in new_coords]) - temporal_param
    t2 = np.max([c[2] for c in new_coords]) + temporal_param
    for i in range(len(full_list)):
        old_event = full_list[i]
        old_coords = old_event[1]
        old_times = [c[2] for c in old_coords]
        time_checks = [t for t in old_times if t >= t1 and t <= t2]

        if len(time_checks) > 0:
            for coord in new_coords:
                # Check if the time coordinate is within an old event
                radii = []
                new_y = coord[0]
                new_x = coord[1]
                for oc in old_coords:
                    old_y = oc[0]
                    old_x = oc[1]
                    dy = abs(old_y - new_y)
                    dx = abs(old_x - new_x)
                    r = np.sqrt((dy ** 2) + (dx ** 2))
                    radii.append(r)
                check = [r for r in radii if r <= radius]
                if any(check):
                    return i, True
                else:
                    return i, False
            else:
                return i, False


def mode(lst):
    return max(set(list(lst)), key=list(lst).count)


def pquery(p, lc, lc_array):
    """Find the landcover code for a particular point (p)."""
    row, col = lc.index(p.x, p.y)
    lc_value = lc_array[row, col]
    return lc_value


def rasterize(src, dst, attribute, resolution, crs, extent, all_touch=False, na=-9999):
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
    cols = int((xmax - xmin)/resolution)
    rows = int((ymax - ymin)/resolution) + 1
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


def requestIO(url):
    """Function for setting IO request for data download"""

    b = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, url)
    c.setopt(c.WRITEFUNCTION, b.write)
    c.perform()
    c.close()
    content = b.getvalue()

    return content


def spCheck(diffs, sp_buf):
    """Quick function to check if events land within the spatial window."""

    checks = [e for e in diffs if abs(e) < sp_buf]
    if any(checks):
        check = True
    else:
        check = False
    return check


def toAcres(p, res):
    return (p*res**2) * 0.000247105


def toHa(p, res):
    return (p*res**2) * 0.0001


def toKms(p, res):
    return (p*res**2)/1000000


def toDays(date, base):
    """Convert dates to days since a base date"""

    if type(date) is str:
        date = dt.datetime.strptime(date, "%Y-%m-%d")
        delta = (date - base)
        days = delta.days
    return days


def asMultiPolygon(polygon):
    if type(polygon) == Polygon:
        polygon = MultiPolygon([polygon])
    return polygon


def get(self, remotepath, localpath=None):
    """Copies a file between the remote host and the local host."""
    """For Paramiko SSH/SFTP Client Access: fuoco.geog.umd.edu"""
    if not localpath:
        localpath = os.path.split(remotepath)[1]
    self._sftp_connect()
    self._sftp.get(remotepath, localpath)


def downloadLC(query, session):
    """Downloads MODIS land cover data"""

    link = query[0]
    dst = query[1]

    try:
        # submit the request using the session
        response = session.get(link, stream=True)
        # raise an exception in case of http errors
        response.raise_for_status()
        # save the file
        with open(dst, 'wb') as fd:
            fd.write(response.content)
    except requests.exceptions.HTTPError as e:
        # handle any errors here
        print(e)


# Classes
class DataGetter:
    """
    Things to do/remember:
        - parallel downloads
    """

    def __init__(self, proj_dir):
        self.proj_dir = proj_dir
        self.date = dt.datetime.today().strftime("%m-%d-%Y")
        self.createPaths()
        self.cpus = os.cpu_count()
        self.modis_template_path = os.path.join(proj_dir, "rasters/")
        self.modis_template_file_root = "mosaic_template.tif"
        self.landcover_path = os.path.join(proj_dir, "rasters/landcover")
        self.landcover_file_root = "lc_mosaic_"
        self.modis_crs = ("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
        self.nc_path = os.path.join(proj_dir, "rasters/burn_area/netcdfs")
        self.hdf_path = os.path.join(proj_dir, "rasters/burn_area/hdfs")
        self.tiles = ["h08v04", "h09v04", "h10v04", "h11v04", "h12v04",
                      "h13v04", "h08v05", "h09v05", "h10v05", "h11v05",
                      "h12v05", "h08v06", "h09v06", "h10v06", "h11v06"]
        print("Project Folder: " + proj_dir)

    def createPaths(self):
        sub_folders = ["rasters/burn_area", "rasters/burn_area/hdfs",
                       "rasters/ecoregion", "rasters/landcover",
                       "rasters/landcover/mosaics/", "shapefiles/ecoregion",
                       "tables"]
        folders = [os.path.join(self.proj_dir, sf) for sf in sub_folders]
        for f in folders:
            if not os.path.exists(f):
                os.makedirs(f)

    def getBurns(self):
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
        ##################################################################################

        # Check into the UMD SFTP fuoco server using Paramiko
        ssh_client = paramiko.SSHClient()
        ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh_client.connect(hostname="fuoco.geog.umd.edu", username="fire", password="burnt")
        print("Connected to 'fuoco.geog.umd.edu' ...")

        # Open the connection to the SFTP
        sftp_client = ssh_client.open_sftp()

        ##################################################################################

        # Use specified tiles or...download all tiles if the list is empty
        if self.tiles[0].lower() != "all":
            tiles = self.tiles
        else:
            sftp_client.chdir('/data/MODIS/C6/MCD64A1/HDF')
            dirs = sftp_client.listdir()
            tiles = dirs

        # print(tiles)

        # Download the available files and catch failed downloads
        for tile in tiles:
            # Find remote folder for the tile
            sftp_folder = '/data/MODIS/C6/MCD64A1/HDF/' + tile

            # Check if remote folder exists and if not, continue to next tile
            try:
                # Change directory to remote folder
                sftp_client.chdir(sftp_folder)

                hdfs = sftp_client.listdir()
                # hdfs = [h for h in hdfs if ".hdf" in h]

                # Make sure local target folder exists
                folder = os.path.join(self.hdf_path, tile)
                if not os.path.exists(folder):
                    os.mkdir(folder)

                # Skip this if the final product exists
                nc_file = os.path.join(self.proj_dir,
                                       "rasters/burn_area/netcdfs/" + tile + ".nc")

                # ~~~~~~~~~~~~~~~~~~Download~~~~~~~~~~~~~~~~~~~~~~~~~~

                if not os.path.exists(nc_file):
                    print("Downloading/Checking HDF files for: " + tile)

                    # Attempt file download
                    try:
                        for h in tqdm(hdfs):
                            remote = sftp_folder+"/"+h
                            os.chdir(folder)
                            sftp_client.get(remote, h)
                    except Exception as e:
                        print(e)

            except Exception:
                print("No MCD64A1 Product for tile: "+str(tile)+", skipping...")

            # Check downloads for missing files
            missings = []
            for hdf in hdfs:
                trgt = os.path.join(folder, hdf)
                remote = os.path.join(sftp_folder, hdf)
                if not os.path.exists(trgt):
                    missings.append(remote)
                else:
                    try:
                        gdal.Open(trgt).GetSubDatasets()[0][0]
                    except Exception:
                        print("Bad file detected, removing to try again...")
                        missings.append(remote)
                        os.remove(trgt)

            # Now try again for the missed files
            if len(missings) > 0:
                print("Missed Files: " + str(missings))
                print("trying again...")

                # # Check into SFTP server again
                # # Should not have to do this twice ...
                # ssh_client = paramiko.SSHClient()
                # ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
                # ssh_client.connect(hostname="fuoco.geog.umd.edu", username="fire", password="burnt")
                #
                # # Open the connection to the SFTP
                # sftp_client = ssh_client.open_sftp()

                for m in missings:
                    tile = m.split("/")[-2]
                    sftp_folder = "/MCD64A1/C6/HDF/" + tile
                    sftp_client.chdir(sftp_folder)
                    file = os.path.basename(m)
                    localpath = os.path.join(self.hdf_path, tile)
                    trgt = os.path.join(self.hdf_path, tile, file)

                    # Attempt re-download
                    try:
                        for h in tqdm(hdfs):
                            remote = sftp_folder+"/"+h
                            os.chdir(localpath)
                            sftp_client.get(remote, localpath)
                    except Exception as e:
                        print(e)

                    # Check the downloaded file
                    try:
                        gdal.Open(trgt).GetSubDatasets()[0][0]
                        missings.remove(file)
                    except Exception as e:
                        print(e)

                # Close new SFTP connection
                sftp_client.close()
                ssh_client.close()

                # If that doesn"t get them all, give up.
                if len(missings) > 0:
                    print("There are still " + str(len(missings)) + " missed files.")
                    print("Try downloading these files manually: ")
                    for mm in missings:
                        print(mm)

        # End connection to Clients
        ssh_client.close()
        sftp_client.close()
        print("Disconnected from 'fuoco.geog.umd.edu' ...")

        # Build the netcdfs here
        tile_files = {}
        for tid in tiles:
            files = glob(os.path.join(self.hdf_path, tid, "*hdf"))
            tile_files[tid] = files

        # Merge one year into a reference mosaic
        if not os.path.exists(self.modis_template_path):
            print("Creating reference mosaic ...")
            folders = glob(os.path.join(self.hdf_path, "*"))
            file_groups = [glob(os.path.join(f, "*hdf")) for f in folders]
            for f in file_groups:
                f.sort()
            files = [f[0] for f in file_groups]
            dss = [rasterio.open(f).subdatasets[0] for f in files]
            tiles = [rasterio.open(d) for d in dss]
            mosaic, transform = merge(tiles)
            crs = tiles[0].meta.copy()
            template_path = os.path.join(self.modis_template_path,
                                         self.modis_template_file_root)
            crs.update({"driver": "GTIFF",
                        "height": mosaic.shape[1],
                        "width": mosaic.shape[2],
                        "transform": transform})

            with rasterio.open(template_path, "w+", **crs) as dst:
                dst.write(mosaic)

        # Build one netcdf per tile
        for tid in tiles:
            files = tile_files[tid]
            if len(files) > 0:
                try:
                    self.buildNCs(files)
                except Exception as e:
                    file_name = os.path.join(self.nc_path, tid + ".nc")
                    print("Error on tile " + tid + ": " + str(e))
                    print("Removing " + file_name + " and moving on.")
                    os.remove(file_name)

    def getEcoregion(self, ecoregion_level=1, rasterize=False):
        # Omernick's Ecoregions - EPA North American Albers
        """
        Update 02/2021: EPA FTP site is glitchy, adding local file in "ref"
        If download fails, use local file
        """
        if not os.path.exists(os.path.join(self.proj_dir, "shapefiles/ecoregion/us_eco_l4.shp")):

            eco_ftp = 'ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip'

            try:
                eco = gpd.read_file(eco_ftp)
                eco.crs = {"init": "epsg:5070"}
                eco.to_file(os.path.join(self.proj_dir,
                                         "shapefiles/ecoregion/NA_CEC_Eco_Level3.shp"))
            except Exception:
                print("Failed to connect to EPA ftp site: using local file ...")
                eco = gpd.read_file(os.path.join(os.getcwd(),
                                    "ref", "us_eco", "NA_CEC_Eco_Level3.shp"))

            ref_cols = ['NA_L3CODE', 'NA_L3NAME',
                        'NA_L2CODE', 'NA_L2NAME',
                        'NA_L1CODE', 'NA_L1NAME',
                        'NA_L3KEY', 'NA_L2KEY', 'NA_L1KEY']

            # Create a reference table for ecoregions
            eco_ref = eco[ref_cols].drop_duplicates()

            # Character cases are inconsistent between I,II and III,IV levels
            def ecoCap(string):
                strings = string.split()
                strings = [s.lower() if s != "USA" else s for s in strings]
                caps = [s.capitalize() if s not in ["USA", "and"] else
                        s for s in strings]
                for i in range(len(caps)):
                    cp = caps[i]
                    if "/" in cp:
                        cp = "/".join([c.capitalize() for c in cp.split("/")])
                    if "-" in cp:
                        cp = "-".join([c.title() for c in cp.split("-")])
                    caps[i] = cp
                return " ".join(caps)

            eco_ref = eco_ref.applymap(ecoCap)
            eco_ref.to_csv(os.path.join(self.proj_dir, "tables/eco_refs.csv"),
                           index=False)

        # Rasterize Omernick Ecoregions
        if rasterize and not os.path.exists(os.path.join(self.proj_dir,
                                            "rasters/ecoregion/NA_CEC_Eco_Level3_modis.tif")):

            # We need something with the correct geometry
            src = eco
            dst = os.path.join(self.proj_dir,
                               "rasters/ecoregion/NA_CEC_Eco_Level3_modis.tif")
            extent_template_file = os.path.join(
                    self.proj_dir, "shapefiles/modis_sinusoidal_grid_world.shp")

            # Getting the extent regardless of existing files from other runs
            template1 = gpd.read_file(extent_template_file)
            template1["h"] = template1["h"].apply(lambda x: "{:02d}".format(x))
            template1["v"] = template1["v"].apply(lambda x: "{:02d}".format(x))
            template1["tile"] = "h" + template1["h"] + "v" + template1["v"]
            template1 = template1[template1["tile"].isin(self.tiles)]

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

            # If these aren"t present, I say just go ahead and download
            exts = []
            for et in extent_tiles:
                folder = os.path.join(self.proj_dir, "rasters/burn_area/hdfs",
                                      et)
                if not os.path.exists(folder):
                    self.getBurns()

                file = glob(os.path.join(folder, "*hdf"))[0]
                file_pointer = gdal.Open(file)
                dataset_pointer = file_pointer.GetSubDatasets()[0][0]
                ds = gdal.Open(dataset_pointer)
                geom = ds.GetGeoTransform()
                ulx, xres, xskew, uly, yskew, yres = geom
                lrx = ulx + (ds.RasterXSize * xres)
                lry = uly + (ds.RasterYSize * yres) + yres
                exts.append([ulx, lry, lrx, uly])

            extent = [exts[0][0], exts[1][1], exts[2][2], exts[3][3]]
            wkt = ds.GetProjection()
            attribute = "US_L3CODE"
            rasterize(src, dst, attribute, xres, wkt, extent)

    def getLandcover(self, landcover_type=1):
        """
        A method to download and process landcover data from NASA"s Land
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

        class SessionWithHeaderRedirection(requests.Session):
            AUTH_HOST = 'urs.earthdata.nasa.gov'

            def __init__(self, username, password):
                super().__init__()
                self.auth = (username, password)

            # Overrides from the library to keep headers
            # When redirected to or from the NASA auth host
            def rebuild_auth(self, prepared_request, response):
                headers = prepared_request.headers
                url = prepared_request.url
                if 'Authorization' in headers:
                    original_parsed = requests.utils.urlparse(response.request.url)
                    redirect_parsed = requests.utils.urlparse(url)
                    if (original_parsed.hostname != redirect_parsed.hostname) and \
                            redirect_parsed.hostname != self.AUTH_HOST and \
                            original_parsed.hostname != self.AUTH_HOST:
                        del headers['Authorization']
                    return

        # Use specified tiles or...
        if self.tiles[0].lower() != "all":
            tiles = self.tiles
        # ...download all tiles if the list is empty
        else:
            # Check into the UMD SFTP fuoco server using Paramiko
            ssh_client = paramiko.SSHClient()
            ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh_client.connect(hostname="fuoco.geog.umd.edu",
                               username="fire", password="burnt")
            print("Connected to 'fuoco.geog.umd.edu' ...")
            # Open the connection to the SFTP
            sftp_client = ssh_client.open_sftp()
            sftp_client.chdir('/data/MODIS/C6/MCD64A1/HDF')
            tiles = sftp_client.listdir()
            ssh_client.close()
            sftp_client.close()

        # Check if the mosaic already exists
        lc_dirs = os.listdir(os.path.join(self.proj_dir, "rasters/landcover/mosaics"))
        if len(lc_dirs) == 0:

            # Get the full string for land cover type
            lp = self.landcover_path
            lc_type = "type" + str(landcover_type)

            # Get available years
            r = requestIO("https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/")
            soup = BeautifulSoup(r, 'html.parser')
            links = [link["href"] for link in soup.find_all("a", href=True)]
            years = [ll[:4] for ll in links if '01.01' in ll]

            # To skip downloaded tiles start by finding local files
            local_files = {os.path.split(fldr)[-1]: files for
                           fldr, _, files in os.walk(self.landcover_path)}
            local_files = {k: v for k, v in local_files.items() if
                           k not in ["landcover", "mosaics"]}

            # Let's check that each file is good, defining function here
            def fileCheck(landcover_path, year, file):
                path = os.path.join(landcover_path, year, file)
                try:
                    rasterio.open(path, 'w+', driver='HDF4Image')
                    return True
                except Exception:
                    return False

            # Check that gdal/rasterio will open the files, if not drop
            for year in local_files.keys():
                files = [f if fileCheck(lp, year, f) else
                         os.remove(os.path.join(lp, year, f))
                         for f in local_files[year]]
                files = [f for f in files if f]
                local_files[year] = files

            # Now, for each local year, get the tiles needed
            needed_tiles = {}
            for year, files in local_files.items():
                local_tiles = [f.split(".")[2] for f in files]
                missing_tiles = [t for t in tiles if t not in local_tiles]
                needed_tiles[year] = missing_tiles

            # Now (this may be empty) fill in all tiles for totally missing years
            for year in years:
                if year not in needed_tiles.keys():
                    needed_tiles[year] = tiles

            # Now, count all of the needed files for each key (year)
            file_count = sum([len(v) for k, v in needed_tiles.items()])

            # If we need anything at all we'll have to do gain access
            if file_count > 0:
                print("Retrieving land cover rasters from NASA's Earthdata " + "service...")
                print("Register at the link below to obtain a username and " + "password:")
                print("https://urs.earthdata.nasa.gov/")
                # Create remote session with user credentials/authentication
                username = input("Enter NASA Earthdata User Name: ")
                password = getpass("Enter NASA Earthdata Password: ")
                session = SessionWithHeaderRedirection(username, password)

                # Get all the remote and local file paths
                queries = []
                print("Retrieving needed landcover file paths...")
                needed_years = [y for y in needed_tiles.keys() if
                                len(needed_tiles[y]) > 0]
                for yr in tqdm(needed_years, position=0, file=sys.stdout):
                    # Make sure destination folder exists
                    if not os.path.exists(os.path.join(lp, yr)):
                        os.mkdir(os.path.join(lp, yr))

                    # Get needed tiles for this year
                    year_tiles = needed_tiles[yr]

                    # Retrieve list of links to hdf files
                    url = ("https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/" + yr + ".01.01/")
                    r = requestIO(url)
                    soup = BeautifulSoup(r, 'html.parser')
                    names = [link["href"] for link in soup.find_all("a", href=True)]
                    names = [n for n in names if "hdf" in n and "xml" not in n]
                    names = [n for n in names if n.split('.')[2] in year_tiles]
                    links = [url + lc for lc in names]

                    # Build list of local file paths and check if they're needed
                    dsts = [os.path.join(lp, yr, names[i]) for
                            i in range(len(links))]

                    # Group links and local paths for parallel downloads
                    query = [(links[i], dsts[i]) for i in range(len(links))
                             if not os.path.exists(dsts[i])]
                    queries = queries + query

                # Now get land cover data from earthdata.nasa.gov
                if len(queries) > 0:
                    print("Retrieving landcover data...")
                    # filename = url[url.rfind('/')+1:]
                    ncores = cpu_count()
                    pool = Pool(int(ncores/2))
                    try:
                        for p in tqdm(pool.map(downloadLC, queries),
                                      total=len(queries), position=0,
                                      file=sys.stdout):
                            p
                    except Exception:
                        try:
                            _ = [downloadLC(q, session) for q in tqdm(queries, position=0, file=sys.stdout)]
                        except Exception as e:
                            template = "Download failed: error type {0}:\n{1!r}"
                            message = template.format(type(e).__name__, e.args)
                            print(message)

            # Now process these tiles into yearly geotiffs. Do this everytime.
            if not os.path.exists(os.path.join(self.proj_dir,
                                               "rasters/landcover/mosaics")):
                os.mkdir(os.path.join(self.proj_dir, "rasters/landcover/mosaics"))

            print("Mosaicking/remosaicking landcover tiles...")
            for y in tqdm(years, position=0, file=sys.stdout):

                # Filter available files for the requested tiles
                lc_files = glob(os.path.join(self.landcover_path, y, "*.hdf"))
                lc_files = [f for f in lc_files if f.split(".")[2] in self.tiles]

                # Use the subdataset name to get the right land cover type
                data_sets = []
                for f in lc_files:
                    subdss = rasterio.open(f).subdatasets
                    trgt_ds = [sd for sd in subdss if lc_type in sd.lower()][0]
                    data_sets.append(trgt_ds)

                    # with rasterio.open(f, driver='HDF4Image') as dataset:
                    #     for name in dataset.subdatasets:
                    #         if lc_type in name.lower():
                    #             data_sets.append(dataset)

                trgt_ds = [sd for sd in subdss if lc_type in sd.lower()][0]
                data_sets.append(trgt_ds)

                # Create pointers to the chosen land cover type
                tiles = [rasterio.open(ds) for ds in data_sets]

                # Mosaic them together
                mosaic, transform = merge(tiles)

                # Get coordinate reference information
                crs = tiles[0].meta.copy()
                crs.update({"driver": "GTIFF",
                            "height": mosaic.shape[1],
                            "width": mosaic.shape[2],
                            "transform": transform})

                # Save mosaic file
                file = self.landcover_file_root + lc_type + "_" + y + ".tif"
                path = os.path.join(self.proj_dir,
                                    "rasters/landcover/mosaics", file)
                with rasterio.open(path, "w+", **crs) as dst:
                    dst.write(mosaic)

            # Print location
            print("Landcover data saved to " + os.path.join(self.proj_dir,
                                                            "rasters/landcover/mosaics"))

    def getShapes(self, ecoregion_level=None):
        """
        Just to grab some basic shapefiles needed for calculating statistics.
        """
        if not os.path.exists(os.path.join(self.proj_dir, "shapefiles")):
            os.mkdir(os.path.join(self.proj_dir, "shapefiles"))

        # Variables
        conus_states = ["WV", "FL", "IL", "MN", "MD", "RI", "ID", "NH", "NC",
                        "VT", "CT", "DE", "NM", "CA", "NJ", "WI", "OR", "NE",
                        "PA", "WA", "LA", "GA", "AL", "UT", "OH", "TX", "CO",
                        "SC", "OK", "TN", "WY", "ND", "KY", "VI", "ME", "NY",
                        "NV", "MI", "AR", "MS", "MO", "MT", "KS", "IN", "SD",
                        "MA", "VA", "DC", "IA", "AZ"]
        modis_crs = self.modis_crs

        # MODIS Sinusoial World Grid
        if not os.path.exists(os.path.join(os.getcwd(), "ref", "modis_sinusoidal_grid_world.shp")):
            print("Downloading MODIS Sinusoidal Projection Grid...")
            src = ("http://book.ecosens.org/wp-content/uploads/2016/06/" + "modis_grid.zip")
            modis = gpd.read_file(src)
            modis.crs = modis_crs
            modis.to_file(os.path.join(self.proj_dir,
                                       "shapefiles/modis_sinusoidal_grid_world.shp"))

        # Contiguous United States - WGS84
        if not os.path.exists(os.path.join(self.proj_dir,
                                           "shapefiles/conus.shp")):
            print("Downloading US state shapefile from the Census Bureau...")
            usa = gpd.read_file("http://www2.census.gov/geo/tiger/GENZ2016/" + "shp/cb_2016_us_state_20m.zip")
            conus = usa[usa["STUSPS"].isin(conus_states)]
            conus.crs = {"init": "epsg:4326", "no_defs": True}
            conus.to_file(os.path.join(self.proj_dir, "shapefiles/conus.shp"))

        # Contiguous United States - MODIS Sinusoidal
        if not os.path.exists(os.path.join(self.proj_dir,
                                           "shapefiles/conus_modis.shp")):
            print("Reprojecting state shapefile to MODIS Sinusoidal...")
            conus = gpd.read_file(os.path.join(self.proj_dir,
                                               "shapefiles/conus.shp"))
            modis_conus = conus.to_crs(modis_crs)
            modis_conus.to_file(os.path.join(self.proj_dir,
                                             "shapefiles/conus_modis.shp"))

    def shapeToTiles(self, shp_path):
        """
        Set or reset the tile list using a shapefile. Where shapes intersect
        with the modis sinusoidal grid determines which tiles to use.
        """

        # MODIS CRS retrieved from a single HDF file
        outCRS = '''PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",
        DATUM["Not specified (based on custom spheroid)",
        SPHEROID["Custom spheroid",6371007.181,0]],
        PRIMEM["Greenwich",0],
        UNIT["degree",0.0174532925199433,
        AUTHORITY["EPSG","9122"]]],
        PROJECTION["Sinusoidal"],
        PARAMETER["longitude_of_center",0],
        PARAMETER["false_easting",0],
        PARAMETER["false_northing",0],
        UNIT["Meter",1],
        AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''

        # Attempt to read in the modis grid and download it if not available
        try:
            grid_path = os.path.join(os.getcwd(), "ref", "modis_grid.gpkg")
            modis_grid = gpd.read_file(grid_path)
            modis_grid.set_crs(outCRS, inplace=True, allow_override=True)
        except Exception:
            print("MODIS Grid not found, downloading from EcoSens...")
            # MODIS Sinusoial World Grid
            src = ("http://book.ecosens.org/wp-content/uploads/2016/06/" + "modis_grid.zip")
            modis = gpd.read_file(src)
            modis_grid.set_crs(outCRS, inplace=True, allow_override=True)
            modis.to_file(os.path.join(self.proj_dir,
                                       "shapefiles", "modis_sinusoidal_grid_world.gpkg"))

        # Read in the input shapefile and reproject to MODIS sinusoidal
        if str(os.path.basename(shp_path).endswith(".shp")):
            try:
                source = gpd.read_file(shp_path)
                source.to_crs(outCRS, inplace=True)
            except Exception as e:
                print("Error: " + str(e))
                print("Failed to reproject file, ensure a coordinate reference system is specified.")
        elif str(os.path.basename(shp_path).endswith(".gpkg")):
            try:
                source = gpd.read_file(shp_path)
                source.to_crs(outCRS, inplace=True)
                print(source.crs)
            except Exception as e:
                print("Error: " + str(e))
                print("Failed to reproject file, ensure a coordinate reference system is specified.")

        # Left join shapefiles with source shape as the left
        shared = gpd.sjoin(source, modis_grid, how="left").dropna()
        shared["h"] = shared["h"].apply(lambda x: "h{:02d}".format(int(x)))
        shared["v"] = shared["v"].apply(lambda x: "v{:02d}".format(int(x)))
        shared["tile"] = shared["h"] + shared["v"]
        tiles = pd.unique(shared["tile"].values)
        self.tiles = tiles

    def buildNCs(self, files):
        """
        Take in a time series of files for the MODIS burn detection dataset and
        create a singular netcdf file.
        """
        savepath = self.nc_path
        # Check that the target folder exists, agian.
        if not os.path.exists(savepath):
            os.mkdir(savepath)

        # Set file names
        files.sort()
        names = [os.path.split(f)[-1] for f in files]
        names = [f.split(".")[2] + "_" + f.split(".")[1][1:] for f in names]
        tile_id = names[0].split("_")[0]
        file_name = os.path.join(savepath, tile_id + ".nc")

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
            y = nco.createVariable("y",  np.float64, ("y",))
            x = nco.createVariable("x",  np.float64, ("x",))
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
                array = convertDates(array, year)
                try:
                    variable[tidx, :, :] = array
                except Exception:
                    print(f + ": failed, probably had wrong dimensions, " + "inserting a blank array in its place.")
                    blank = np.zeros((ny, nx))
                    variable[tidx, :, :] = blank
                tidx += 1

            # Done
            nco.close()


class EventPerimeter:
    def __init__(self, event_id, coord_list=[]):
        self.event_id = event_id
        self.merge_id = np.nan
        self.coords = []
        self.coords = self.add_coordinates(coord_list)

    def add_coordinates(self, coord_list):
        for coord in coord_list:
            self.coords.append(coord)
        return self.coords

    def get_event_id(self):
        return self.event_id

    def get_merge_id(self):
        return self.merge_id

    def get_coords(self):
        return self.coords


class EventGrid:
    """
    For a single file, however large, find sites with any burn detections
    in the study period, loop through these sites and group by the space-time
    window, save grouped events to a data frame on disk.
    """

    def __init__(self, proj_dir, nc_path=("rasters/burn_area/netcdfs/"),
                 spatial_param=5, temporal_param=11, area_unit="Unknown",
                 time_unit="days since 1970-01-01"):
        self.proj_dir = proj_dir
        self.nc_path = os.path.join(proj_dir, nc_path)
        self.spatial_param = spatial_param
        self.temporal_param = temporal_param
        self.area_unit = area_unit
        self.time_unit = time_unit
        self.event_grid = {}
        self.next_event_id = 1
        self.input_array = self.get_input_xarray()

    def get_input_xarray(self):
        burns = xr.open_dataset(self.nc_path)
        self.data_set = burns
        self.coordinates = burns.coords
        input_array = burns.value

        return input_array

    def add_event_grid(self, event_id, new_pts):
        for p in new_pts:
            entry = {p: event_id}
            self.event_grid.update(entry)

    def merge_perimeters(self, perimeters, event_id, obsolete_id):
        # set the merge id in the obsolete id
        perimeters[obsolete_id-1].merge_id = event_id
        new_pts = perimeters[obsolete_id-1].coords

        # update the event_grid and add points to event_id perimeter
        for p in new_pts:
            self.event_grid[p] = event_id
        perimeters[event_id-1].add_coordinates(new_pts)

        # set old perimeter to null
        merge_notice = "Merged with event {}".format(event_id)
        perimeters[obsolete_id-1].coords = [merge_notice, new_pts]

        return perimeters

    def get_spatial_window(self, y, x, array_dims):
        """
        Pull in the spatial window around a detected event and determine its
        shape and the position of the original point within it. Finding this
        origin point is related to another time saving step in the event
        classification procedure.
        """
        top = max(0, y - self.spatial_param)
        bottom = min(array_dims[0], y + self.spatial_param)
        left = max(0, x - self.spatial_param)
        right = min(array_dims[1], x + self.spatial_param)

        #  Derive full xarray coordinates from just the window for speed
        ydim = array_dims[0]
        xdim = array_dims[1]

        # Expand the spatial dimension
        tps = [i for i in range(self.spatial_param)]

        # There are four edge cases
        x_edges_1 = [0 + t for t in tps]
        x_edges_2 = [xdim - t for t in tps]
        y_edges_1 = [0 + t for t in tps]
        y_edges_2 = [ydim - t for t in tps]

        # Get the full y, x coords of the origin, and window coords of date
        if y in y_edges_1:
            ycenter = y
            oy = 0
        elif y in y_edges_2:
            ycenter = y - ydim
            oy = y - self.spatial_param
        else:
            ycenter = self.spatial_param
            oy = y - self.spatial_param
        if x in x_edges_1:
            xcenter = x
            ox = 0
        elif x in x_edges_2:
            xcenter = x - xdim
            ox = x - self.spatial_param
        else:
            xcenter = self.spatial_param
            ox = x - self.spatial_param
        center = [ycenter, xcenter]
        origin = [oy, ox]

        return top, bottom, left, right, center, origin

    def get_availables(self):
        """
        To save time, avoid checking cells with no events at any time step.
        Create a mask of max values at each point. If the maximum at a cell is
        less than or equal to zero there were no values and it will not be
        checked in the event classification step.
        """
        # Low memory - Somehow leads to slow loop in get_event_perimeters
        # We want to get the mask without pulling the whole thing into memory
#        burns = xr.open_dataset(self.nc_path, chunks={"x": 500, "y": 500})
#
#        # Pull in only the single max value array
#        mask = burns.max(dim="time").compute()
#
#        # Get the y, x positions where one or more burns were detected
#        locs = np.where(mask.value.values > 0)
#
#        # Now pair these
#        available_pairs = []
#        for i in range(len(locs[0])):
#            available_pairs.append([locs[0][i], locs[1][i]])
#
#        # Leaving the data set open causes problems
#        burns.close()

        # Using memory - can handle large tiles, but gets pretty high
        mask = self.input_array.max(dim="time")
        locs = np.where(mask > 0)
        available_pairs = []
        for i in range(len(locs[0])):
            available_pairs.append([locs[0][i], locs[1][i]])
        del mask

        # Done.
        return available_pairs

    def get_event_perimeters(self):
        """
        Iterate through each cell in the 3D MODIS Burn Date tile and group it
        into fire events using the space-time window.
        """
        print("Filtering out cells with no events...")
        available_pairs = self.get_availables()

        # If the other pointers aren"t closed first, this will be very slow
        arr = self.input_array

        # This is to check the window positions
        nz, ny, nx = arr.shape
        dims = [ny, nx]
        perimeters = []

        # traverse spatially, processing each burn day
        print("Building event perimeters...")
        for pair in tqdm(available_pairs, position=0, file=sys.stdout):
            # Separate coordinates
            y, x = pair

            # get the spatial window
            [top, bottom, left,
             right, center, origin] = self.get_spatial_window(y, x, dims)
            cy, cx = center

            # what if we pull in the window?
            window = arr[:, top:bottom+1, left:right+1].data

            # The center of the window is the target burn day
            center_burn = window[:, cy, cx]
            center_burn = center_burn[center_burn > 0]

            # Loop through each event in the window and identify neighbors
            for burn in center_burn:
                new_pts = []
                curr_event_ids = []

                # Now we can get the values and position right away
                diff = abs(burn - window)
                val_locs = np.where(diff <= self.temporal_param)
                y_locs = val_locs[1]
                x_locs = val_locs[2]
                oy, ox = origin

                # Get the actual x, y positions from the window coordinates
                vals = window[val_locs]
                ys = [oy + yl for yl in y_locs]
                xs = [ox + xl for xl in x_locs]

                # Now get the geographic coordinates from tile positions
                all_ys = self.coordinates["y"].data
                all_xs = self.coordinates["x"].data
                ys = all_ys[ys]
                xs = all_xs[xs]

                # Now check if this point is in the event_grid yet
                for i in range(len(vals)):
                    curr_pt = (float(ys[i]), float(xs[i]), float(vals[i]))

                    # already assigned to an event
                    if (curr_pt in self.event_grid):
                        if self.event_grid[curr_pt] not in curr_event_ids:
                            curr_event_ids.append(self.event_grid[curr_pt])
                    else:
                        new_pts.append(curr_pt)

                # If this is a new event
                if len(curr_event_ids) == 0:
                    # create a new perimeter object
                    perimeter = EventPerimeter(self.next_event_id, new_pts)

                    # append to perimeters list
                    perimeters.append(perimeter)

                    # add points to the grid
                    self.add_event_grid(self.next_event_id, new_pts)

                    # increment the event ID
                    self.next_event_id += 1

                # If all points part of same existing event
                elif len(curr_event_ids) == 1:
                    event_id = curr_event_ids[0]
                    if len(new_pts):
                        perimeters[event_id - 1].add_coordinates(new_pts)
                        self.add_event_grid(event_id, new_pts)

                # events overlap
                else:
                    perimeters = self.merge_perimeters(perimeters,
                                                       curr_event_ids[0],
                                                       curr_event_ids[1])

        return perimeters


class ModelBuilder:
    def __init__(self, file_name, proj_dir, tiles, daily, shapefile, spatial_param=5,
                 temporal_param=11, landcover_type=None,
                 ecoregion_type=None, ecoregion_level=None):
        self.file_name = file_name
        self.proj_dir = proj_dir
        self.tiles = tiles
        self.spatial_param = spatial_param
        self.temporal_param = temporal_param
        self.landcover_type = landcover_type
        self.ecoregion_type = ecoregion_type
        self.ecoregion_level = ecoregion_level
        self.daily = daily
        self.shapefile = shapefile
        self.getFiles(file_name)
        self.setGeometry()

    def getFiles(self, file_name):
        # Get the requested files names
        files = []
        for t in self.tiles:
            path = os.path.join(
                    self.proj_dir, "rasters/burn_area/netcdfs/" + t + ".nc")
            files.append(path)
        files.sort()
        self.files = files

        # Make sure the main data frame file name has the right extension
        if ".csv" not in file_name:
            self.file_name = os.path.splitext(file_name)[0] + ".csv"
        else:
            self.file_name = self.file_name

    def setGeometry(self):
        # Use the first file to extract some more parameters for later
        builder = EventGrid(nc_path=self.files[0],
                            proj_dir=self.proj_dir,
                            spatial_param=self.spatial_param,
                            temporal_param=self.temporal_param)
        self.crs = builder.data_set.crs
        self.geom = self.crs.geo_transform
        self.res = self.geom[1]
        self.sp_buf = builder.spatial_param * self.res

    def buildEvents(self):
        """
        Use the EventGrid class to classify events tile by tile and then merge
        them all together for a seamless set of wildfire events.
        """
        # Make sure the destination folder exists
        if not(os.path.exists(os.path.dirname(self.file_name))):
            print(self.file_name)
            os.makedirs(os.path.dirname(self.file_name))

        # Create an empty list and data frame for the
        tile_list = []
        columns = ["id", "date", "x", "y", "edge", "tile"]
        df = pd.DataFrame(columns=columns)
        base = dt.datetime(1970, 1, 1)

        # Loop through each netcdf file and build individual tile events
        for file in self.files:
            tile_id = file[-9:-3]
            if os.path.exists(
                    os.path.join(
                        self.proj_dir, "tables/events/" + tile_id + ".csv")
                    ):
                print(tile_id + " event table exists, skipping...")
            elif not os.path.exists(
                    os.path.join(self.proj_dir, "rasters/burn_area/netcdfs/" + tile_id + ".nc")):
                pass
            else:
                print("\n" + tile_id)

                # Create a new event object
                builder = EventGrid(nc_path=file,
                                    proj_dir=self.proj_dir,
                                    spatial_param=self.spatial_param,
                                    temporal_param=self.temporal_param)

                # Classify event perimeters
                perims = builder.get_event_perimeters()

                # Remove empty perimeters
                perims = [p for p in perims if type(p.coords[0]) is not str]
                tile_list.append(perims)

                # Extract just the event ID, days, and x,y MODIS coordinates
                plist = [(p.get_event_id(), p.coords) for p in perims]

                # Identify edge cases, so either x or y is within 5 cells
                maxys = builder.data_set["y"].data[:builder.spatial_param]
                minys = builder.data_set["y"].data[-builder.spatial_param:]
                maxxs = builder.data_set["x"].data[-builder.spatial_param:]
                minxs = builder.data_set["x"].data[:builder.spatial_param]
                yedges = list(maxys) + list(minys)
                xedges = list(maxxs) + list(minxs)

                # Create an empty data frame
                print("Building data frame...")
                events = []
                coords = []
                edges = []
                ys = []
                xs = []
                dates = []
                for p in plist:
                    coord = [list(c) for c in p[1]]
                    edge = [edgeCheck(yedges, xedges, c, self.sp_buf) for
                            c in coord]
                    if any(edge):
                        edge = [True for e in edge]
                    event = list(np.repeat(p[0], len(coord)))
                    y = [c[0] for c in coord]
                    x = [c[1] for c in coord]
                    date = [base + dt.timedelta(c[2]) for c in coord]
                    events.append(event)
                    coords.append(coord)
                    edges.append(edge)
                    ys.append(y)
                    xs.append(x)
                    dates.append(date)

                # Flatten each list of lists
                events = flttn(events)
                coords = flttn(coords)
                edges = flttn(edges)
                ys = flttn(ys)
                xs = flttn(xs)
                dates = flttn(dates)
                edf = pd.DataFrame(
                        OrderedDict({"id": events, "date": dates, "x": xs,
                                     "y": ys, "edge": edges, "tile": tile_id})
                        )
                if not os.path.exists(os.path.join(self.proj_dir, "tables/events")):
                    os.mkdir(os.path.join(self.proj_dir, "tables/events"))
                edf.to_csv(
                    os.path.join(
                        self.proj_dir, "tables/events/" + tile_id + ".csv"),
                    index=False)

        # Clear memory
        gc.collect()

        # Now read in the event data frames (use dask, instead, save memory)
        print("Reading saved event tables back into memory...")
        efiles = glob(os.path.join(self.proj_dir, "tables/events/*csv"))
        efiles = [e for e in efiles if e[-10:-4] in self.tiles]
        edfs = [pd.read_csv(e) for e in efiles]

        # Merge with existing records
        print("Concatenating event tables...")
        df = pd.concat(edfs)

        def toDays(date, base):
            date = dt.datetime.strptime(date, "%Y-%m-%d")
            delta = (date - base)
            days = delta.days
            return days

        print("Creating unique ids...")
        df["id"] = df["tile"] + "_" + df["id"].astype(str)

        print("Converting days since 1970 to dates...")
        df["days"] = df["date"].apply(toDays, base=base)

        # Cut the edge events out into a separate df
        print("Separating tile edge events from center events...")
        edges = df[df["edge"] == True]
        not_edges = df[df["edge"] == False]

        # Merge where needed
        print("Merging edge-case tile events...")
        eids = list(edges["id"].unique())
        for iden in tqdm(eids, position=0, file=sys.stdout):
            # Split, one vs all
            edf = edges[edges["id"] == iden]
            edf2 = edges[edges["id"] != iden]
            days = edf["days"]

            # Sometimes these are empty
            try:
                d1 = min(days)
                d2 = max(days)
            except Exception:
                pass

            # If events aren't close enough in time the list will be empty
            edf2 = edf2[(abs(edf2["days"] - d1) < self.temporal_param) | (abs(edf2["days"] - d2) < self.temporal_param)]
            eids2 = list(edf2["id"].unique())

            # If there are event close in time, are they close in space?
            for iden2 in eids2:
                edf2 = edges[edges["id"] == iden2]
                ydiffs = [y - edf2["y"].values for y in edf["y"].values]
                xdiffs = [x - edf2["x"].values for x in edf["x"].values]
                ychecks = [spCheck(yds, self.sp_buf) for yds in ydiffs]
                xchecks = [spCheck(xds, self.sp_buf) for xds in xdiffs]
                checks = [ychecks[i] * xchecks[i] for i in range(len(ychecks))]
                if any(checks):
                    # Merge events! Merge into the earliest event
                    # MC - Need to replace this with the event builder on these points (?)
                    d12 = edf2["days"].min()
                    if d1 < d12:
                        edges["id"][edges["id"] == iden2] = iden
                    else:
                        edges["id"][edges["id"] == iden] = iden2

        # Concatenate edge df back into main df
        print("Recombining edge and center cases...")
        df = pd.concat([not_edges, edges])

        # Reset id values in chronological order
        print("Resetting ids in chronological order..")
        df["first"] = df.groupby("id").days.transform("min")
        firsts = df[["id", "first"]].drop_duplicates()
        firsts = firsts.sort_values("first")
        firsts["new_id"] = range(1, firsts.shape[0] + 1)
        idmap = dict(zip(firsts["id"], firsts["new_id"]))
        df["id"] = df["id"].map(idmap)
        df = df.sort_values("id")

        # put these in order
        df = df[["id", "tile", "date", "x", "y"]]

        # Finally save
        print("Saving merged event table to " + self.file_name)
        df.to_csv(self.file_name, index=False)

    def buildPoints(self):
        '''
        Build point GeoDataFrame from classified fire event table
        '''
        # Read in the event table
        print("Reading classified fire event table...")
        df = pd.read_csv(self.file_name, low_memory=False)

        # Go ahead and create daily id (did) for later
        df["did"] = df["id"].astype(str) + "-" + df["date"].astype(str)

        # Get geometries
        crs = self.crs
        geom = self.geom
        proj4 = crs.proj4
        resolutions = [geom[1], geom[-1]]

        # Filter columns, center pixel coordinates, and remove repeating pixels
        df["x"] = df["x"] + (resolutions[0]/2)
        df["y"] = df["y"] + (resolutions[1]/2)

        # Each entry gets a point object from the x and y coordinates.
        print("Converting data frame to spatial object...")
        df["geometry"] = df[["x", "y"]].apply(lambda x: Point(tuple(x)),
                                              axis=1)

        gdf = gpd.GeoDataFrame(df, crs=proj4, geometry=df["geometry"])

        return gdf

    def buildFireAttributes(self):
        '''
        Take the data table, add in attributes, and overwrite file.
        '''
        # If there is an event table first, make a spatial object out of it
        if os.path.exists(self.file_name):
            gdf = self.buildPoints()
        else:
            print("Run BuildEvents first...")
            return

        # Space is tight and we need the spatial resolution
        res = self.res

        # Adding fire attributes
        print("Adding fire attributes ...")

        # gdf["did"] = gdf["id"].astype(str) + "-" + gdf["date"].astype(str)
        gdf['pixels'] = gdf.groupby(['id', 'date'])['id'].transform('count')

        group = gdf.groupby('id')

        # gdf['max_growth_dates'] = group[['date', 'pixels']].apply(maxGrowthDate)
        gdf['date'] = gdf['date'].apply(lambda x: dt.datetime.strptime(x, '%Y-%m-%d'))

        gdf['ig_date'] = group['date'].transform('min')
<<<<<<< HEAD
        gdf['ig_day'] = gdf['ig_date'].apply(lambda x: dt.datetime.strftime(x, '%j'))
=======
        gdf['ig_day'] = gdf['ign_date'].apply(lambda x: dt.datetime.strftime(x, '%j'))
>>>>>>> f9c937947a7242699bce38cc69e415131d30db38
        gdf['ig_month'] = gdf['ig_date'].apply(lambda x: x.month)
        gdf['ig_year'] = gdf['ig_date'].apply(lambda x: x.year)
        gdf['last_date'] = group['date'].transform('max')

        gdf['tot_pix'] = group['id'].transform('count')

        gdf['daily_duration'] = gdf['date'] - gdf['ig_date']
        gdf['event_day'] = gdf['daily_duration'].apply(lambda x: x.days + 1)
        gdf['final_duration'] = gdf['last_date'] - gdf['ig_date']
        gdf['event_dur'] = gdf['final_duration'].apply(lambda x: x.days + 1)
        gdf.drop('final_duration', axis=1)

        gdf['dy_ar_km2'] = gdf['pixels'].apply(toKms, res=res)
        gdf['tot_ar_km2'] = gdf['tot_pix'].apply(toKms, res=res)

        # # Add in cumulative sum attributes
        # gdf['cml_pixels'] = gdf.groupby('id')['pixels'].transform(pd.Series.cumsum)
        # gdf['cml_area_km2'] = gdf['cml_pixels'].apply(toKms, res=res)
        # gdf['perc_total_area_km2'] = (gdf['daily_area_km2'] / gdf['total_area_km2'] * 100).astype(int)
        # gdf['perc_cml_area_km2'] = (gdf['cml_area_km2'] / gdf['total_area_km2'] * 100).astype(int)

        gdf['fsr_px_dy'] = gdf['tot_pix'] / gdf['event_dur']
        gdf['fsr_km2_dy'] = gdf['fsr_px_dy'].apply(toKms, res=res)

        gdf['mx_grw_px'] = group['pixels'].transform('max')
        gdf['mn_grw_px'] = group['pixels'].transform('min')
        gdf['mu_grw_px'] = group['pixels'].transform('mean')

        gdf['mx_grw_km2'] = gdf['mx_grw_px'].apply(toKms, res=res)
        gdf['mn_grw_km2'] = gdf['mn_grw_px'].apply(toKms, res=res)
        gdf['mu_grw_km2'] = gdf['mu_grw_px'].apply(toKms, res=res)

        max_date = pd.DataFrame(group[['date', 'pixels']].apply(maxGrowthDate).reset_index())
        max_date = max_date.rename(columns={0: 'mx_grw_dte'})
        gdf = pd.merge(gdf, max_date, on="id")

        gdf = gdf[['id', 'date', 'ig_date', 'ig_day', 'ig_month',
                   'ig_year', 'last_date', 'event_day', 'event_dur',
                   'pixels', 'tot_pix', 'dy_ar_km2', 'tot_ar_km2',
                   'fsr_px_dy', 'fsr_km2_dy',
                   'mx_grw_px', 'mn_grw_px', 'mu_grw_px',
                   'mx_grw_km2', 'mn_grw_km2', 'mu_grw_km2', 'mx_grw_dte',
                   'x', 'y', 'geometry']]

        gdf = gdf.reset_index(drop=True)

        ############################################
        # Retrieve landcover attributes if requested
        if self.landcover_type is not None:
            print('Adding landcover attributes...')

            # We'll need to specify which type of landcover
            lc_types = {1: "IGBP global vegetation classification scheme",
                        2: "University of Maryland (UMD) scheme",
                        3: "MODIS-derived LAI/fPAR scheme",
                        4: "MODIS-derived Net Primary Production (NPP) scheme",
                        5: "Plant Functional Type (PFT) scheme."}

            # Get mosaicked landcover geotiffs
            lc_files = glob(os.path.join(self.proj_dir,
                                         "rasters/landcover/mosaics",
                                         "*tif"))
            lc_files.sort()
            lc_years = [int(f[-8:-4]) for f in lc_files]
            lc_files = {lc_years[i]: lc_files[i] for i in range(len(lc_files))}

            # Rasterio point querier (will only work here)
            def pointQuery(row):
                x = row['x']
                y = row['y']
                try:
                    val = int([val for val in lc.sample([(x, y)])][0])
                except:
                    val = np.nan
                return val

            # Get the range of burn years
            burn_years = list(gdf['ig_year'].unique())

            # This works faster when split by year and the pointer is outside
            # This is also not the best way
            sgdfs = []
            for year in tqdm(burn_years, position=0,
                             file=sys.stdout):

                sgdf = gdf[gdf['ig_year'] == year]
<<<<<<< HEAD

                # Now set year one back for landcover
                year = year-1
=======
>>>>>>> f9c937947a7242699bce38cc69e415131d30db38

                # Use previous year's lc
                if year < min(lc_years):
                    year = min(lc_years)
                    # print("Min: ", year)
                elif year > max(lc_years):
                    year = max(lc_years)
                    # print("Max: ", year)

                lc_file = lc_files[year]
                lc = rasterio.open(lc_file)
                sgdf['lc_code'] = sgdf.apply(pointQuery, axis=1)
                sgdfs.append(sgdf)

            gdf = pd.concat(sgdfs)
            gdf = gdf.reset_index(drop=True)
            gdf['lc_mode'] = gdf.groupby('id')['lc_code'].transform(mode)
            # gdf = gdf.drop("landcover", axis=1)
            # Add in the class description from landcover tables
            lc_table = pd.read_csv(os.path.join(self.proj_dir, 'tables', 'landcover',
                                                'MCD12Q1_LegendDesc_Type{}.csv'.format(str(self.landcover_type))))
            gdf = pd.merge(left=gdf, right=lc_table, how='left', left_on='lc_mode', right_on='Value')
            gdf = gdf.drop('Value', axis=1)
            gdf['lc_type'] = lc_types[int(self.landcover_type)]
<<<<<<< HEAD
            gdf.rename({'lc_description': 'lc_desc'}, inplace=True, axis='columns')
=======
>>>>>>> f9c937947a7242699bce38cc69e415131d30db38

        ############################################
        # Retrieve ecoregion attributes if requested
        # Add the ecoregion attributes
        if self.ecoregion_type or self.ecoregion_level:
            print("Adding ecoregion attributes ...")
            if self.ecoregion_level or self.ecoregion_type == "na":
                # Different levels have different sources
                eco_types = {
                    'NA_L3CODE': ('Level III Ecoregions ' + '(NA-Commission for Environmental Cooperation)'),
                    'NA_L2CODE': ('Level II Ecoregions ' + '(NA-Commission for Environmental Cooperation)'),
                    'NA_L1CODE': ('Level I Ecoregions ' + '(NA-Commission for Environmental Cooperation)')}

                # Read in the Level File (contains every level) and reference table
                shp_path = os.path.join(self.proj_dir, "shapefiles/ecoregion/NA_CEC_Eco_Level3.gpkg")
                eco = gpd.read_file(shp_path)
                eco.to_crs(gdf.crs, inplace=True)

                # Filter for selected level (level III defaults to US-EPA version)
                if not self.ecoregion_level:
                    self.ecoregion_level = 3

                eco_code = [c for c in eco.columns if str(self.ecoregion_level) in
                            c and 'CODE' in c]

                if len(eco_code) > 1:
                    eco_code = [c for c in eco_code if 'NA' in c][0]
                else:
                    eco_code = eco_code[0]

                print("Selected ecoregion code: "+str(eco_code))

                # Find modal eco region for each event id
                # print(eco[eco.columns[-10:]].head())
                eco = eco[[eco_code, 'geometry']]
                gdf = gpd.sjoin(gdf, eco, how="left", op="within")
                gdf = gdf.reset_index(drop=True)
                gdf['eco_mode'] = gdf.groupby('id')[eco_code].transform(mode)

                # gdf[eco_code] = gdf[eco_code].apply(
                #        lambda x: int(x) if not pd.isna(x) else np.nan)

                # Add in the type of ecoregion
                gdf['eco_type'] = eco_types[eco_code]

                # Add in the name of the modal ecoregion
                eco_ref = pd.read_csv(os.path.join(self.proj_dir,
                                                   'tables/eco_refs.csv'))
                eco_name = eco_code.replace('CODE', 'NAME')
                eco_df = eco_ref[[eco_code, eco_name]].drop_duplicates()
                eco_map = dict(zip(eco_df[eco_code], eco_df[eco_name]))
                gdf['eco_name'] = gdf[eco_code].map(eco_map)

                # Clean up column names
                gdf = gdf.drop('index_right', axis=1)
                gdf = gdf.drop(eco_code, axis=1)
                # gdfd.rename({eco_code: 'eco_mode'}, inplace=True, axis='columns')

            else:

                # Read in the world ecoregions from WWF
                eco_path = os.path.join(os.getcwd(), "ref", "world_ecoregions", "wwf_terr_ecos_modis.shp")
                eco = gpd.read_file(eco_path)
                eco = eco.to_crs(gdf.crs, inplace=True)

                # Find modal eco region for each event id
                gdf = gpd.sjoin(gdf, eco, how="left", op="within")
                gdf = gdf.reset_index(drop=True)

                gdf["eco_mode"] = gdf.groupby('id')['ECO_NUM'].transform(mode)
                gdf["eco_name"] = gdf["ECO_NAME"]
                gdf["eco_type"] = "WWF Terrestrial Ecoregions of the World"

        # Save event level attributes
        print("Overwriting data frame at " + self.file_name + "...")
        gdf.to_csv(self.file_name, index=False)

    def buildPolygons(self, daily_shp_path, event_shp_path):

        # Make sure we have the target folders
        if not(os.path.exists(os.path.dirname(event_shp_path))):
            os.makedirs(os.path.dirname(event_shp_path))

        # grab the modis crs and space
        res = self.res

        # Create a spatial points object
        gdf = self.buildPoints()
        # gdf.to_crs(crs, inplace=True)

        # Create a circle buffer
        print("Creating buffer...")
        geometry = gdf.buffer(1 + (res/2))
        gdf["geometry"] = geometry

        # Then create a square envelope around the circle
        gdf["geometry"] = gdf.envelope

        # Now add the first date of each event and merge daily event detections
        print("Dissolving polygons...")
        gdfd = gdf.dissolve(by="did", as_index=False)

        # Add the ignition coords, convert to lat/long
        gdfd_x = pd.DataFrame(gdfd.groupby('id')['x'].nth(0))
        gdfd_y = pd.DataFrame(gdfd.groupby('id')['y'].nth(0))
        gdfd_xy = gdfd_x.merge(gdfd_y, on='id')
        gdfd = gdfd.merge(gdfd_xy, on='id')
        gdfd = gdfd.rename(columns={"x_x": "x", "y_x": "y",
                                    "x_y": "ig_utm_x",
                                    "y_y": "ig_utm_y"})
        gdfd = gdfd.drop(['x', 'y'], axis=1)

        # Cast as multipolygon for each geometry
        print("Converting polygons to multipolygons...")
        gdfd["geometry"] = gdfd["geometry"].apply(asMultiPolygon)

        # Save the daily before dissolving into event level
        # Only save the daily polygons if user specified to do so
        if self.daily == "yes":
            print("Saving daily file to " + daily_shp_path)
            gdfd.to_csv(str(self.file_name)[:-4]+"_daily"+".csv", index=False)
            if self.shapefile:
                # gdf.to_crs(outCRS, inplace=True)
                gdfd.to_file(daily_shp_path, driver="GPKG")

        # Drop the daily attributes before exporting event-level
        gdf = gdfd.drop(['did', 'pixels', 'date', 'event_day',
                        'dy_ar_km2'], axis=1)

        # Dissolve by ID to create event-level
        gdf = gdf.dissolve(by="id", as_index=False)

        # Calculate perimeter length
        print("Calculating perimeter lengths...")
        gdf["tot_perim"] = gdf["geometry"].length

        # We still can't have multiple polygon types
        print("Converting polygons to multipolygons...")
        gdf["geometry"] = gdf["geometry"].apply(asMultiPolygon)

        # Export event-level to CSV
        gdf.to_csv(str(self.file_name)[:-4]+"_events"+".csv", index=False)
        # Save as gpkg if specified
        if self.shapefile:
            # Now save as a geopackage and csv
            print("Saving event-level file to " + event_shp_path)
            # gdf.to_crs(outCRS, inplace=True)
            gdf.to_file(event_shp_path, driver="GPKG")

        # Remove the intermediate file
        os.remove(self.file_name)
