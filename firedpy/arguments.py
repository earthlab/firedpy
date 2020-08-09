#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Just a script for sample arguments.

Created on Tue Oct 22 08:37:01 2019

@author: travis
"""

import argparse
from functions import DataGetter, ModelBuilder
import os
import tempfile
import time
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# Start the timer (seconds)
start = time.perf_counter()

# Call help statements
data_help = ("""
    The project directory you would like to use for input and output
    data files. Defaults to a temporary directory.
    """)
file_help = ("""
    The file name of the resulting dataframe. This will be saved in
    the "outputs/tables" folder of the chosen project directory. Defaults
    to "modis_events.csv".
    """)
eco_help = ("""
    Provide this option to associate each event with an ecoregion_level.
    """)
lc_help = ("""
        Provide this option if you would like to associate each event with
        a land cover category. If so, you will have to register at NASA's
        Earthdata service (https://urs.earthdata.nasa.gov/home) and enter your
        user name and password when prompted. Land cover comes from the
        MODIS/Terra+Aqua Land Cover Type data set (MCD12Q1). Enter a number for
        one of the following land cover types:
            1: IGBP global vegetation classification scheme,
            2: University of Maryland (UMD) scheme,
            3: MODIS-derived LAI/fPAR scheme,
            4: MODIS-derived Net Primary Production (NPP) scheme,
            5: Plant Functional Type (PFT) scheme.
        (e.g., -landcover 4). Defaults to none.
        """)
shp_help = ("""
    Provide this option if you would like to build shapefiles from the
    event data frame. Shapefiles of both daily progression and overall
    event perimeters will be written to the "outputs/shapefiles" folder of
    the chosen project directory. These will be saved in geopackage format
    (.gpkg) using the file basename of the fire event data frame (e.g.
    'modis_events_daily.gpkg' and 'modis_events.gpkg')
    """)
sp_help = ("""
    The number of cells (default MODIS resolution is 463 meters) to search for
    neighboring burn detections. Defaults to 5 cells.
    """)
tile_help = ("""
    You may specify the tiles as a list of characters (no quotes no spaces)
    (e.g., h08v04 h09v04 ...) or leave this blank to default to tiles
    covering the Contiguous United States. Specify "all" to use all
    available MODIS tiles. Alternatively, provide a path to a shapefile
    with either a ".shp" or ".gpkg" extension to use intersecting MODIS
    tiles.
    """)
tmp_help = ("""
    The number of days to search for neighboring burn detections. Defaults
    to 11 days between events.
    """)

# Provide arguments
parser = argparse.ArgumentParser()
parser.add_argument("-proj_dir", dest="proj_dir",
                    default=tempfile.mkdtemp(), help=data_help)
parser.add_argument("-file_name", dest="file_name",
                    default="modis_events.csv",
                    help=file_help)
parser.add_argument("-ecoregion_level", dest="ecoregion_level", default=4,
                    help=lc_help)
parser.add_argument("-landcover_type", dest="landcover_type", default=1,
                    help=lc_help)
parser.add_argument("--shapefile", action='store_true', help=shp_help)
parser.add_argument("-spatial_param", dest="spatial_param", default=5,
                    type=int, help=sp_help)
parser.add_argument("-temporal_param", dest="temporal_param", default=11,
                    type=int, help=tmp_help)
parser.add_argument("-tiles", "--names-list", nargs="+", dest="tiles",
                    default=["h08v04"],
                        help=tile_help)

# Parse argument responses
args = parser.parse_args()
proj_dir = '/home/travis/fired' #args.proj_dir
file_name = os.path.join(proj_dir, "outputs", "tables", args.file_name)
ecoregion_level = args.ecoregion_level
landcover_type = args.landcover_type
spatial_param = args.spatial_param
temporal_param = args.temporal_param
tiles = args.tiles
shapefile = args.shapefile

# Make sure the project directory exists
if not os.path.exists(proj_dir):
    os.makedirs(proj_dir)

# Create data object
data = DataGetter(proj_dir)

# Assign target MODIS tiles to the data object
if os.path.splitext(tiles[0])[1] in [".shp", ".gpkg"]:
    shp = tiles[0]
    print("Filtering for MODIS tiles that intersect \n    " + shp)
    data.shapeToTiles(shp)
    tiles = data.tiles
else:
    data.tiles = tiles

# Create Model Builder object
models = ModelBuilder(file_name=file_name,
                      proj_dir=proj_dir,
                      tiles=tiles,
                      spatial_param=spatial_param,
                      temporal_param=temporal_param,
                      landcover_type=landcover_type,
                      ecoregion_level=ecoregion_level)

# Use arguments for shapefile source and destination file paths
file_base = os.path.splitext(os.path.basename(file_name))[0]
daily_shp_file = "_".join([file_base, "daily"])
daily_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                              daily_shp_file + ".gpkg")
event_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                              file_base + ".gpkg")
