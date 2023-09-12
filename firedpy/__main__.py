# -*- coding: utf-8 -*-
import argparse
import os
#import shutil
import time
import warnings
#import sys
from typing import List, Any

#from http.cookiejar import CookieJar
#import urllib.request
#from .functions import DataGetter, ModelBuilder
#from utilities.create_readme import makeReadMe

warnings.filterwarnings("ignore", category=FutureWarning)

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))

# Call help statements
DATA_HELP = ("""
    The project directory you would like to use for  and output
    data files. Defaults to a temporary directory 'firedpy/proj'.
    """)
FILE_HELP = ("""
    The file name of the resulting dataframe. This will be saved in
    the "outputs/tables" folder of the chosen project directory. Defaults
    to "fired_events.csv" and "fired_daily.csv" if daily data is requested.
    """)
DAILY_HELP = ("""
    You may specify whether to create the daily polygons or just the event-level perimeter
    for your analysis area. Options are "yes" (to create the daily polygons and the event polygons),
    "no" (create the event level only).
    """)
ECO_HELP = ("""
    You can specify the ecoregion type as either "world" or "na":
    "world" = World Terrestrial Ecoregions (World Wildlife Fund (WWF))
    "na" = North American ecoregions (Omernick, 1987)

    Most common (modal) ecoregion across the event is used.

    Further, to associate each event with North American ecoregions (Omernick,
    1987) you may provide a number corresponding to an ecoregion level. Ecoregions
    are retrieved from www.epa.gov and levels I through IV are available.
    Levels I and II were developed by the North American Commission for
    Environmental Cooperation. Levels III and IV were developed by the
    United States Environmental Protection Agency. For events with more
    than one ecoregion, the most common value will be used. Defaults to
    none.
    """)
LC_HELP = ("""
    To include land cover as an attribute, provide a number corresponding
    with a MODIS/Terra+Aqua Land Cover (MCD12Q1) category followed with :username:password of your NASA's Earthdata service
    account. Available land cover categories:
        1: IGBP global vegetation classification scheme,
        2: University of Maryland (UMD) scheme,
        3: MODIS-derived LAI/fPAR scheme.

    If you do not have an account register at https://urs.earthdata.nasa.gov/home. Defaults to none.
    """)
SHP_HELP = ("""
    Provide this option if you would like to build shapefiles from the
    event data frame. Spefify either "shp", "gpkg", or both. Shapefiles of both daily progression and overall
    event perimeters will be written to the "outputs/shapefiles" folder of
    the chosen project directory. These will be saved in the specified geopackage format
    (.gpkg), ERSI Shapefile format (.shp), or save them in both formats using the file basename of the fire event data frame (e.g.
    'modis_events_daily.gpkg' and 'modis_events.gpkg')
    """)
SP_HELP = ("""
    The number of cells (~463 m resolution) to search for neighboring burn
    detections. Defaults to 5 cells in all directions.
    """)
TILE_HELP = ("""
    You may specify the tiles as a list of characters (no quotes no spaces)
    (e.g., h08v04 h09v04 ...) or leave this blank to default to tiles
    covering the Contiguous United States. Specify "all" to use all
    available MODIS tiles. Alternatively, provide a path to a shapefile
    with either a ".shp" or ".gpkg" extension to use intersecting MODIS
    tiles. In the firedpy directory, you can access any of the 50 states by specifying ref/us_states/state_name.gpkg,
    all 7 continents by specifying ref/continents/continent_name.gpkg, and a country by ref/individual_countries/country_name.gpkg.
    A list of all avalible countries can be found in the ReadMe.
    """)
TMP_HELP = ("""
    The number of days to search for neighboring burn detections. Defaults
    to 11 days between events.
    """)
START_YR = ("""
    The first year of fired events.
    """)
END_YR = ("""
    The last year of fired events.
    """)
FULL_CSV = ("""
    If included full attribute table will exported to csv. If not included only x and y coordinates, event date, and
    event id will be exported to a csv.""")


class FiredpyArgumentParser:

    def __init__(self, params_file):
        self.params_file = params_file
        self.arguments = self._load_params()

    def _load_params(self):
        args = {}
        try:
            with open(self.params_file, 'r') as file:
                for line in file:
                    name, prompt, arg_type, last_value, accepted_values = line.strip().split(',')
                    args[name] = {
                        'prompt': prompt.replace('\\n', '\n'),
                        'type': arg_type,
                        'last_value': last_value if last_value != 'none' else None,
                        'accepted_values': accepted_values if accepted_values != 'none' else accepted_values.split('|')
                    }
        except FileNotFoundError:
            pass

        return args

    def _save_params(self):
        with open(self.params_file, 'w') as file:
            for name, data in self.arguments.items():
                file.write(f"{name},f{data['prompt']},{data['type']},{data['last_value']},{data['accepted_values']}\n")
        self._load_params()

    def prompt_for_argument(self, arg_name, prompt_override: str = None, accepted_value_override: List[Any] = None):
        if arg_name not in self.arguments:
            raise ValueError(f"Argument {arg_name} not found in params file.")

        arg_data = self.arguments[arg_name]
        last_value = arg_data['last_value']
        prompt = arg_data['prompt'] if prompt_override is None else prompt_override
        accepted_values = arg_data['accepted_values'] if accepted_value_override is None else accepted_value_override
        user_input = input(f"{prompt} [Default: {last_value}]: ")

        # If no input is given, use the last value.
        if not user_input:
            user_input = last_value

        if accepted_values is not None and user_input not in accepted_values:
            print(f'{user_input} is not an acceptable value. Acceptable values are:\n {accepted_values}')
            return self.prompt_for_argument(arg_name, prompt_override)

        # Convert the input to the appropriate type.
        if arg_data['type'] == 'int':
            user_input = int(user_input)
        elif arg_data['type'] == 'float':
            user_input = float(user_input)
        elif arg_data['type'] == 'bool':
            user_input = user_input.lower() in ['y', 'true', 'yes']

        # Update the last used value.
        self.arguments[arg_name]['last_value'] = user_input

        # Save updated parameters to file.
        self._save_params()

        return user_input


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

    daily = firedpy_parser.prompt_for_argument('daily')
    spatial_param = firedpy_parser.prompt_for_argument('spatial_param')
    temporal_param = firedpy_parser.prompt_for_argument('temporal_param')
    shp_type = firedpy_parser.prompt_for_argument('shp_type')
    shapefile = shp_type != 'none'
    if not shapefile:
        shp_type = None
    file_name = "fired_" + str(tile_name)
    file_path = os.path.join(out_dir, "outputs", "tables", file_name)

    print("""""")
    landcover_type = input("")
    if landcover_type != '':
        landcover_type = int(landcover_type)
        print(
            "To get the landcover you need to have a username and password with NASA Earthdata services. You can register at the link below to obtain a username and " + "password:")
        print("https://urs.earthdata.nasa.gov/")
        username = input("Enter NASA Earthdata User Name: ")
        password = input("Enter NASA Earthdata Password: ")

    elif landcover_type == '':
        landcover_type = None
        username = ''
        password = ''

    start_yr = input("Enter the year you want to start or press enter for all dates: ")
    end_yr = input("Enter the year you want to end or press enter for all dates: ")
    if start_yr != '' and end_yr != '':
        start_yr = int(start_yr)
        end_yr = int(end_yr)
    else:
        start_yr = None
        end_yr = None
    temp = 1
    full_csv = input("Enter True if you want a full csv. Enter False if you want a raw csv. ")
#     else:
#
#         # Parse argument responses
#         args = parser.parse_args()
#         proj_dir = args.proj_dir
#         ecoregion_type = args.ecoregion_type
#         ecoregion_level = args.ecoregion_level
#         landcover_type = args.landcover_type
#         daily = args.daily
#         spatial_param = args.spatial_param
#         temporal_param = args.temporal_param
#         tiles = args.tiles
#         shp_type = args.shp_type
#         start_yr = args.start_yr
#         end_yr = args.end_yr
#         full_csv = args.full_csv
#         if shp_type != 'none':
#             shapefile = True
#         elif shp_type == 'none':
#             shapefile = False
#             shp_type = None
#         if landcover_type:
#             user_pass = landcover_type.split(':', 2)
#             username = str(user_pass[1])
#             password = str(user_pass[2])
#             landcover_type = int(user_pass[0])
#         else:
#             username = ''
#             password = ''
#         if tiles:
#             temp = 2
#             name = str(tiles[0])
#             nums = [str(0), str(1), str(4), str(3), str(4), str(5), str(6), str(7), str(8), str(9)]
#             for i in nums:
#                 if i in name:
#                     temp = 3
#                     tilename = tiles
#             if temp != 3:
#                 check_file = tiles[0]
#                 while True:
#                     if os.path.exists(check_file):
#                         break
#                     else:
#                         print("Not a valid file name.")
#                         check_file = input("Enter aoi again in ref/dir_name/aoi_file.gpkg format:")
#                         tiles = [check_file, ]
#                         name = str(tiles[0])
#                 name = name.split("/")
#                 name = name[-1]
#                 name = name.split(".")
#                 tilename = name[0]
#                 file_path = os.path.join(args.proj_dir,
#                                          "outputs", "tables",
#                                          args.file_name + "_" + str(tilename))
#             else:
#                 # Assign the temporary file name including the spatial and temporal parameters
#                 file_path = os.path.join(args.proj_dir,
#                                          "outputs", "tables",
#                                          args.file_name)
#         else:
#             temp = 2
#             tilename = "CONUS"
#             file_path = os.path.join(args.proj_dir,
#                                      "outputs", "tables",
#                                      args.file_name)
#
#             # Assign the temporary file name including the spatial and temporal parameters
#
#     # Transfer the lookup tables
#     if landcover_type:
#         # Earthdata Login
#         # test url for correct user/password
#         url = "https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2019.01.01/MCD12Q1.A2019001.h13v12.006.2020212130349.hdf"
#
#         password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
#         password_manager.add_password(None, "https://urs.earthdata.nasa.gov", username, password)
#         # Create a cookie jar for storing cookies. This is used to store and return
#         # the session cookie given to use by the data server (otherwise it will just
#         # keep sending us back to Earthdata Login to authenticate).  Ideally, we
#         # should use a file based cookie jar to preserve cookies between runs. This
#         # will make it much more efficient.
#         cookie_jar = CookieJar()
#         # Install all the handlers.
#         opener = urllib.request.build_opener(
#             urllib.request.HTTPBasicAuthHandler(password_manager),
#             # urllib.request.HTTPHandler(debuglevel=1),    # Uncomment these two lines to see
#             # urllib.request.HTTPSHandler(debuglevel=1),   # details of the requests/responses
#             urllib.request.HTTPCookieProcessor(cookie_jar))
#         urllib.request.install_opener(opener)
#         ##Checking to make sure username and password is correct:
#         check = None
#         while check is None:
#             try:
#                 request = urllib.request.Request(url)
#                 response = urllib.request.urlopen(request)
#                 check = 1
#             except Exception:
#                 print("Invalid username or password for NASA Earthdata service account. Try again \n")
#                 # Try again
#                 username = input("Enter NASA Earthdata User Name: ")
#                 password = input("Enter NASA Earthdata Password: ")
#                 password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
#                 password_manager.add_password(None, "https://urs.earthdata.nasa.gov", username, password)
#                 cookie_jar = CookieJar()
#                 opener = urllib.request.build_opener(
#                     urllib.request.HTTPBasicAuthHandler(password_manager),
#                     # urllib.request.HTTPHandler(debuglevel=1),    # Uncomment these two lines to see
#                     # urllib.request.HTTPSHandler(debuglevel=1),   # details of the requests/responses
#                     urllib.request.HTTPCookieProcessor(cookie_jar))
#                 urllib.request.install_opener(opener)
#
#         lookup = os.path.join(os.getcwd(), 'ref', 'landcover',
#                               'MCD12Q1_LegendDesc_Type{}.csv'.format(str(landcover_type)))
#
#         new_path = os.path.join(proj_dir, 'tables', 'landcover')
#         if not os.path.exists(new_path):
#             os.makedirs(new_path)
#         new_file = os.path.join(new_path, 'MCD12Q1_LegendDesc_Type{}.csv'.format(str(landcover_type)))
#         shutil.copy(lookup, new_file)
#     # Make sure the project directory exists
#     if not os.path.exists(proj_dir):
#         os.makedirs(proj_dir)
#
#     # Get ecoregions if requested, use local file first
#     if ecoregion_type or ecoregion_level:
#
#         new_path = os.path.join(proj_dir, 'shapefiles', 'ecoregion')
#
#         if not os.path.exists(new_path):
#             os.makedirs(new_path)
#
#         if ecoregion_type == 'world':
#             fname = 'wwf_terr_ecos.gpkg'
#             lookup = os.path.join(os.getcwd(), 'ref', 'world_ecoregions', fname)
#         elif ecoregion_type == 'na' or ecoregion_level:
#             fname = 'NA_CEC_Eco_Level3.gpkg'
#             lookup = os.path.join(os.getcwd(), 'ref', 'us_eco', fname)
#         try:
#             new_file = os.path.join(new_path, fname)
#             shutil.copy(lookup, new_file)
#         except Exception:
#             # data.getEcoregion(ecoregion_level)
#             pass
#
#     # Create data object
#     data = DataGetter(proj_dir, start_yr, end_yr, username, password)
#
#     data.getEcoregion(ecoregion_level)
#
#     # Assign target MODIS tiles to the data object
#     if os.path.splitext(tiles[0])[1] in [".shp", ".gpkg"]:
#         shp = tiles[0]
#         print("Filtering for MODIS tiles that intersect \n    " + shp)
#         data.shapeToTiles(shp)
#         tiles = data.tiles
#     else:
#         data.tiles = tiles
#         shp = tiles[0]
#
#     # Get all of the MODIS burn area hdfs
#     try:
#         data.getBurns()
#     except (KeyboardInterrupt, SystemExit):
#         raise
#     except Exception as e:
#         template = "\nDownload failed: error type {0}:\n{1!r}"
#         message = template.format(type(e).__name__, e.args)
#         print(message)
#
#     # Download land cover if requested
#     if landcover_type is not None:
#         data.getLandcover(landcover_type)
#
#     # Grab the basename for the output file name
#     file_base = os.path.basename(file_path)
#
#     # Add date range to the file names before exporting final data frame
#     date_range = []
#     # TODO: Should be using regex and file type class here
#     for root, dirs, files in os.walk(os.path.join(proj_dir, 'rasters', 'burn_area', 'hdfs')):
#         for f in files:
#             if "MCD64A1" in f:
#                 dr = int(f.split('.')[1][1:])
#                 date_range.append(dr)
#     last_date = sorted(date_range)[-1]
#     first_date = sorted(date_range)[1]
#     file_name = os.path.join(os.path.dirname(file_path), file_base + "_to" + str(last_date))
#
#     # Create Model Builder object
#     models = ModelBuilder(file_name=file_name,
#                           proj_dir=proj_dir,
#                           tiles=tiles,
#                           shp=shp,
#                           spatial_param=spatial_param,
#                           temporal_param=temporal_param,
#                           landcover_type=landcover_type,
#                           ecoregion_type=ecoregion_type,
#                           ecoregion_level=ecoregion_level,
#                           daily=daily,
#                           shapefile=shapefile,
#                           shp_type=shp_type)
#
#     # Now go ahead and create the events (Memory's a bit tight for parallel)
#     models.buildEvents()
#
#     # Now add fire attributes to this table
#     models.buildFireAttributes()
#
#     # And now build the polygons
#     pref = file_base + "_to" + str(last_date)
#     daily_shp_file = "_".join([pref, "daily"])
#     event_shp_file = "_".join([pref, "events"])
#
#     if shp_type == "shp":
#         daily_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
#                                           daily_shp_file + ".shp")
#         event_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
#                                           event_shp_file + ".shp")
#     if shp_type == "both":
#         daily_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
#                                           daily_shp_file + ".shp")
#         daily_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
#                                       daily_shp_file + ".gpkg")
#         event_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
#                                           event_shp_file + ".shp")
#         event_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
#                                       event_shp_file + ".gpkg")
#     else:
#         daily_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
#                                       daily_shp_file + ".gpkg")
#         event_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
#                                       event_shp_file + ".gpkg")
#         daily_shp_path_shp = ''
#         event_shp_path_shp = ''
#
#     models.buildPolygons(daily_shp_path=daily_shp_path,
#                          event_shp_path=event_shp_path,
#                          daily_shp_path_shp=daily_shp_path_shp,
#                          event_shp_path_shp=event_shp_path_shp,
#                          full_csv=full_csv)
#     makeReadMe(proj_dir, tilename, file_base, temp, first_date, last_date, ecoregion_type, ecoregion_level,
#                landcover_type, daily, spatial_param, temporal_param, shapefile, shp_type)
#     # Print the time it took
#     end = time.perf_counter()
#     seconds = end - start
#     minutes = seconds / 60
#     print("Job completed in {} minutes".format(round(minutes, 2)))
#
#
# if __name__ == "__main__":
#     main()
