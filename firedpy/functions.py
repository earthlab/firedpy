# -*- coding: utf-8 -*-
from bs4 import BeautifulSoup
from collections import OrderedDict
import datetime as dt
import gc
import geopandas as gpd
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
from typing import List, Union
from scipy.spatial import cKDTree


PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))


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

def date_range(perimeter):
    """Converts days in a perimeter object since Jan 1 1970 to date strings"""

    if len(perimeter.spacetime_coordinates) > 0:
        base = dt.datetime(1970, 1, 1)
        days = [p[2] for p in perimeter.spacetime_coordinates]
        day1 = (base + dt.timedelta(days=int(min(days)))).strftime("%Y-%m-%d")
    else:
        day1 = "N/A"
    return day1





def flatten(lst):
    """Just a quick way to flatten lists of lists"""
    lst = [l for sl in lst for l in sl]
    return lst


def max_growth_date(x):
    dates = x["date"].to_numpy()
    pixels = x["pixels"].to_numpy()
    loc = np.where(pixels == np.max(pixels))[0]
    d = np.unique(dates[loc])[0]
    return d


def merge_checker(new_coords, full_list, temporal_param, radius):
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


def request_io(url):
    """Function for setting IO request for data download"""
    b = BytesIO()
    c = pycurl.Curl()
    c.setopt(c.URL, url)
    c.setopt(c.WRITEFUNCTION, b.write)
    c.perform()
    c.close()
    content = b.getvalue()

    return content


def sp_check(diffs, sp_buf):
    """Quick function to check if events land within the spatial window."""

    checks = [e for e in diffs if abs(e) < sp_buf]
    if any(checks):
        check = True
    else:
        check = False
    return check


def to_acres(p, res):
    return (p * res ** 2) * 0.000247105


def to_ha(p, res):
    return (p * res ** 2) * 0.0001


def to_kms(p, res):
    return (p * res ** 2) / 1000000


def to_days(date, base):
    """Convert dates to days since a base date"""

    if type(date) is str:
        date = dt.datetime.strptime(date, "%Y-%m-%d")
        delta = (date - base)
        days = delta.days
    return days


def as_multi_polygon(polygon):
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


