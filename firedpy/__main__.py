# -*- coding: utf-8 -*-
import argparse
import os
#import shutil
import time
import warnings
#import sys
from typing import List, Any

from http.cookiejar import CookieJar
import urllib.request
from .functions import DataGetter, ModelBuilder
from utilities.create_readme import makeReadMe
from firedpy import *
from utilities.argument_parser import FiredpyArgumentParser
from firedpy.spatial import shape_to_tiles


warnings.filterwarnings("ignore", category=FutureWarning)

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))


def get_tile_name_from_directory(parser: FiredpyArgumentParser, directory: str, prompt_message: str):
    """
    Function to prompt the user for a tile name from a given directory.
    """
    available_files = [file.replace('.gpkg', '') for file in os.listdir(directory) if
                       file.endswith('.gpkg')]
    return parser.prompt_for_argument(
        arg_name='tile_name',
        prompt_override=prompt_message,
        accepted_value_override=available_files
    )


def main():
    # Start the timer (seconds)
    start = time.perf_counter()

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-out_dir", dest="out_dir", default=os.path.join(os.getcwd(), 'output'),
                        help=DATA_HELP, required=False)
    parser.add_argument("-file_name", dest="file_name", default="fired", help=FILE_HELP, required=False)
    parser.add_argument("-ecoregion_type", dest="ecoregion_type", help=ECO_HELP, required=False)
    parser.add_argument("-ecoregion_level", dest="ecoregion_level", type=int, default=None, help=ECO_HELP)
    parser.add_argument("-landcover_type", dest="landcover_type", default=None, help=LC_HELP)
    parser.add_argument("-shp_type", dest="shp_type", help=SHP_HELP, default="gpkg")
    parser.add_argument("-spatial", dest="spatial_param", default=5,
                        type=int, help=SP_HELP)
    parser.add_argument("-temporal", dest="temporal_param", default=11,
                        type=int, help=TMP_HELP)
    parser.add_argument("-aoi", "--names-list", nargs="+", dest="tiles",
                        default=["h08v04", "h09v04", "h10v04", "h11v04",
                                 "h12v04", "h13v04", "h08v05", "h09v05",
                                 "h10v05", "h11v05", "h13v04", "h08v05",
                                 "h09v05", "h10v05", "h11v05", "h12v05",
                                 "h08v06", "h09v06", "h10v06", "h11v06"],
                        help=TILE_HELP)
    parser.add_argument("-daily", dest="daily", default="no", help=DAILY_HELP)
    parser.add_argument("-start_yr", dest="start_yr", type=int, default=None, help=START_YR)
    parser.add_argument("-end_yr", dest="end_yr", type=int, default=None, help=END_YR)
    parser.add_argument("--full_csv", action='store_true', help=FULL_CSV)
    args = parser.parse_args()

    firedpy_parser = FiredpyArgumentParser(os.path.join(PROJECT_DIR, 'data', 'params.txt'))

    # Resolve the output directory
    out_dir = args.out_dir if args.out_dir is not None else firedpy_parser.prompt_for_argument('out_dir')
    os.makedirs(out_dir, exist_ok=True)

    # Resolve full csv
    full_csv = args.full_csv if args.full_csv is True else firedpy_parser.prompt_for_argument('full_csv')

    # Resolve tile choice
    tile_choice = firedpy_parser.prompt_for_argument('tile_choice')
    if tile_choice == 'a':
        tile_name = get_tile_name_from_directory(firedpy_parser, os.path.join(PROJECT_DIR, 'ref', 'continents'),
                                                 "Please enter the continent name:")
        tiles = [os.path.join("ref", "continents", tile_name + ".gpkg")]
        ecoregion_type, ecoregion_level = ('na' if tile_name == "north_america" else 'world', 3
        if tile_name == "north_america" else None)

    elif tile_choice == 'b':
        tile_name = get_tile_name_from_directory(firedpy_parser,
                                                 os.path.join(PROJECT_DIR, 'ref', 'individual_countries'),
                                                 "Please enter the country name:")
        tiles = [os.path.join('ref', 'individual_countries', tile_name + ".gpkg")]
        ecoregion_type, ecoregion_level = ('na', 3) if tile_name in [
            'united_States_of_america', 'canada', 'united_states_virgin_islands'] else ('world', None)

    elif tile_choice == 'c':
        tile_name = get_tile_name_from_directory(firedpy_parser, os.path.join(PROJECT_DIR, 'ref', 'us_states'),
                                                 "Please enter the state name:")
        tiles = [os.path.join('ref', 'us_states', tile_name + ".gpkg")]
        ecoregion_type, ecoregion_level = 'na', 3

    elif tile_choice == 'd':
        tiles = firedpy_parser.prompt_for_argument(
            arg_name='tile_name',
            prompt_override="Please enter tiles as a list of characters (no quotes no spaces)(e.g., h08v04 h09v04 "\
                            "...):").split(' ')
        # TODO: Some type of check here for valid tile names?
        tile_name = tiles[0]
        ecoregion_type = firedpy_parser.prompt_for_argument('ecoregion_type')
        ecoregion_level = 3 if ecoregion_type == 'na' else None

    else:
        raise ValueError('Invalid tile choice')

    daily = firedpy_parser.prompt_for_argument('daily')
    spatial_param = firedpy_parser.prompt_for_argument('spatial_param')
    temporal_param = firedpy_parser.prompt_for_argument('temporal_param')

    shp_type = firedpy_parser.prompt_for_argument('shp_type')
    shapefile = shp_type != 'none'
    if not shapefile:
        shp_type = None

    file_name = "fired_" + str(tile_name)
    file_path = os.path.join(out_dir, "outputs", "tables", file_name)

    # Resolve land cover type
    land_cover_type = firedpy_parser.prompt_for_argument('land_cover_type')
    if land_cover_type != 0:
        username = os.environ.get('FIREDPY_ED_USER', None)
        password = os.environ.get('FIREDPY_ED_PWD', None)
        if username is None or password is None:
            print("Please input your NASA Earthdata username and password in order to download the land cover data. If"
                  " you do not have an Earthdata account, you can register at https://urs.earthdata.nasa.gov/. To "
                  "avoid seeing this prompt again, you can set the FIREDPY_ED_USER and FIREDPY_ED_PWD environment "
                  " variables.")
            username = firedpy_parser.prompt_for_argument('username')
            password = firedpy_parser.prompt_for_argument('password', sensitive=True)

    else:
        land_cover_type = None
        username = ''
        password = ''

    start_year = firedpy_parser.prompt_for_argument('start_year')
    end_year = firedpy_parser.prompt_for_argument('end_year')
    if start_year == 'all':
        start_year = None
    if end_year == 'all':
        end_year = None

    temp = 1

    # Transfer the lookup tables
    if land_cover_type is not None:
        # Earthdata Login
        # test url for correct user/password
        url = "https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2019.01.01/MCD12Q1.A2019001.h13v12.006.2020212130349.hdf"

        password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
        password_manager.add_password(None, "https://urs.earthdata.nasa.gov", username, password)
        # Create a cookie jar for storing cookies. This is used to store and return
        # the session cookie given to use by the data server (otherwise it will just
        # keep sending us back to Earthdata Login to authenticate).  Ideally, we
        # should use a file based cookie jar to preserve cookies between runs. This
        # will make it much more efficient.
        cookie_jar = CookieJar()
        # Install all the handlers.
        opener = urllib.request.build_opener(
            urllib.request.HTTPBasicAuthHandler(password_manager),
            # urllib.request.HTTPHandler(debuglevel=1),    # Uncomment these two lines to see
            # urllib.request.HTTPSHandler(debuglevel=1),   # details of the requests/responses
            urllib.request.HTTPCookieProcessor(cookie_jar))
        urllib.request.install_opener(opener)
        ##Checking to make sure username and password is correct:
        check = None
        while check is None:
            try:
                request = urllib.request.Request(url)
                response = urllib.request.urlopen(request)
                check = 1
            except Exception:
                print("Invalid username or password for NASA Earthdata service account. Try again \n")
                # Try again
                username = input("Enter NASA Earthdata User Name: ")
                password = input("Enter NASA Earthdata Password: ")
                password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
                password_manager.add_password(None, "https://urs.earthdata.nasa.gov", username, password)
                cookie_jar = CookieJar()
                opener = urllib.request.build_opener(
                    urllib.request.HTTPBasicAuthHandler(password_manager),
                    # urllib.request.HTTPHandler(debuglevel=1),    # Uncomment these two lines to see
                    # urllib.request.HTTPSHandler(debuglevel=1),   # details of the requests/responses
                    urllib.request.HTTPCookieProcessor(cookie_jar))
                urllib.request.install_opener(opener)

        lookup = os.path.join(os.getcwd(), 'ref', 'landcover',
                              'MCD12Q1_LegendDesc_Type{}.csv'.format(str(land_cover_type)))

        new_path = os.path.join(out_dir, 'tables', 'landcover')
        if not os.path.exists(new_path):
            os.makedirs(new_path)
        new_file = os.path.join(new_path, 'MCD12Q1_LegendDesc_Type{}.csv'.format(str(land_cover_type)))
        shutil.copy(lookup, new_file)

    # Get ecoregions if requested, use local file first
    if ecoregion_type or ecoregion_level:

        new_path = os.path.join(out_dir, 'shapefiles', 'ecoregion')

        if not os.path.exists(new_path):
            os.makedirs(new_path)

        if ecoregion_type == 'world':
            fname = 'wwf_terr_ecos.gpkg'
            lookup = os.path.join(os.getcwd(), 'ref', 'world_ecoregions', fname)
        elif ecoregion_type == 'na' or ecoregion_level:
            fname = 'NA_CEC_Eco_Level3.gpkg'
            lookup = os.path.join(os.getcwd(), 'ref', 'us_eco', fname)
        try:
            new_file = os.path.join(new_path, fname)
            shutil.copy(lookup, new_file)
        except Exception:
            # data.getEcoregion(ecoregion_level)
            pass

    # Assign target MODIS tiles to the data object
    shape_file_path = tiles[0]
    if shape_file_path.endswith(('.shp', '.gpkg')):
        print("Filtering for MODIS tiles that intersect \n    " + shape_file_path)
        tiles = shape_to_tiles(shape_file_path)

    data = DataGetter(out_dir, start_year, end_year, username, password, tiles=tiles)
    data.get_eco_region(ecoregion_level)

    # Get all of the MODIS burn area hdfs
    try:
        data.getBurns()
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        template = "\nDownload failed: error type {0}:\n{1!r}"
        message = template.format(type(e).__name__, e.args)
        print(message)

    # Download land cover if requested
    if landcover_type is not None:
        data.getLandcover(landcover_type)

    # Grab the basename for the output file name
    file_base = os.path.basename(file_path)

    # Add date range to the file names before exporting final data frame
    date_range = []
    # TODO: Should be using regex and file type class here
    for root, dirs, files in os.walk(os.path.join(proj_dir, 'rasters', 'burn_area', 'hdfs')):
        for f in files:
            if "MCD64A1" in f:
                dr = int(f.split('.')[1][1:])
                date_range.append(dr)
    last_date = sorted(date_range)[-1]
    first_date = sorted(date_range)[1]
    file_name = os.path.join(os.path.dirname(file_path), file_base + "_to" + str(last_date))

    # Create Model Builder object
    models = ModelBuilder(file_name=file_name,
                          proj_dir=proj_dir,
                          tiles=tiles,
                          shp=shp,
                          spatial_param=spatial_param,
                          temporal_param=temporal_param,
                          landcover_type=landcover_type,
                          ecoregion_type=ecoregion_type,
                          ecoregion_level=ecoregion_level,
                          daily=daily,
                          shapefile=shapefile,
                          shp_type=shp_type)

    # Now go ahead and create the events (Memory's a bit tight for parallel)
    models.buildEvents()

    # Now add fire attributes to this table
    models.buildFireAttributes()

    # And now build the polygons
    pref = file_base + "_to" + str(last_date)
    daily_shp_file = "_".join([pref, "daily"])
    event_shp_file = "_".join([pref, "events"])

    if shp_type == "shp":
        daily_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
                                          daily_shp_file + ".shp")
        event_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
                                          event_shp_file + ".shp")
    if shp_type == "both":
        daily_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
                                          daily_shp_file + ".shp")
        daily_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      daily_shp_file + ".gpkg")
        event_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
                                          event_shp_file + ".shp")
        event_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      event_shp_file + ".gpkg")
    else:
        daily_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      daily_shp_file + ".gpkg")
        event_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      event_shp_file + ".gpkg")
        daily_shp_path_shp = ''
        event_shp_path_shp = ''

    models.buildPolygons(daily_shp_path=daily_shp_path,
                         event_shp_path=event_shp_path,
                         daily_shp_path_shp=daily_shp_path_shp,
                         event_shp_path_shp=event_shp_path_shp,
                         full_csv=full_csv)
    makeReadMe(proj_dir, tilename, file_base, temp, first_date, last_date, ecoregion_type, ecoregion_level,
               landcover_type, daily, spatial_param, temporal_param, shapefile, shp_type)
    # Print the time it took
    end = time.perf_counter()
    seconds = end - start
    minutes = seconds / 60
    print("Job completed in {} minutes".format(round(minutes, 2)))


if __name__ == "__main__":
    main()
