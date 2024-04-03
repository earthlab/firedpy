import argparse
import os
import shutil
import time
import warnings
import resource

from http.cookiejar import CookieJar
import urllib.request
from src.data_classes import EcoRegion, BurnData, LandCover
from src.model_classes import ModelBuilder
from utilities.create_readme import make_read_me
from src import *
from utilities.argument_parser import FiredpyArgumentParser
from src.spatial import shape_to_tiles
from src.enums import EcoRegionType, TileChoice, ShapeType, LandCoverType

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


def str_to_bool(s: str):
    accepted_values = ['true', 'y', 'yes', 'false', 'n', 'no']
    if s.lower() not in accepted_values:
        raise ValueError(f'{s} must be in {accepted_values}')
    return s.lower() in ['true', 'y', 'yes']


def test_earthdata_credentials(username: str, password: str) -> None:
    # Earthdata Login
    # test url for correct user/password
    url = "https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.061/2019.01.01/BROWSE.MCD12Q1.A2019001.h10v09.061.2022169160720.1.jpg"

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

    request = urllib.request.Request(url)
    urllib.request.urlopen(request)


def peak_memory():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


def cleanup_intermediate_files(out_dir):
    shutil.rmtree(os.path.join(out_dir, 'rasters', 'burn_area'))
    shutil.rmtree(os.path.join(out_dir, 'rasters', 'land_cover'))


def main():
    initial_memory = peak_memory()
    # Start the timer (seconds)
    start = time.perf_counter()

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-out_dir", dest="out_dir", default=os.path.join(PROJECT_DIR, 'output'),
                        help=DATA_HELP)
    parser.add_argument("-file_name", dest="file_name", default="fired", type=str, help=FILE_HELP)
    parser.add_argument("-eco_region_type", type=EcoRegionType, help=ECO_HELP)
    parser.add_argument("-eco_region_level", type=int, help=ECO_HELP)
    parser.add_argument("-land_cover_type", help=LC_HELP, type=int)
    parser.add_argument("-shape_type", help=SHP_HELP, type=ShapeType)
    parser.add_argument("-spatial", dest="spatial_param", type=int, help=SP_HELP)
    parser.add_argument("-temporal", dest="temporal_param", type=int, help=TMP_HELP)
    parser.add_argument("-aoi", "--names-list", nargs="+", dest="tiles",
                        default=["h08v04", "h09v04", "h10v04", "h11v04",
                                 "h12v04", "h13v04", "h08v05", "h09v05",
                                 "h10v05", "h11v05", "h13v04", "h08v05",
                                 "h09v05", "h10v05", "h11v05", "h12v05",
                                 "h08v06", "h09v06", "h10v06", "h11v06"],
                        help=TILE_HELP)
    parser.add_argument("--tile_choice", type=TileChoice, help='')
    parser.add_argument('--tile_name', type=str, help='The name of the tile you would like to choose based'
                                                      ' on your tile choice')
    parser.add_argument("--daily", type=str, help=DAILY_HELP)
    parser.add_argument("-start_year", type=int, help=START_YR)
    parser.add_argument("-end_year", type=int, help=END_YR)
    parser.add_argument("--full_csv", type=str, help=FULL_CSV)
    parser.add_argument('--n_cores', type=int, help='Number of cores to use for parallel processing.')
    parser.add_argument('--cleanup', type=str,
                        help='If set then the burn area and landcover files will be removed after each run to save ' \
                             ' disk space in between multiple runs')
    args = parser.parse_args()

    firedpy_parser = FiredpyArgumentParser(os.path.join(PROJECT_DIR, 'data', 'params.txt'))

    # Resolve the output directory
    out_dir = args.out_dir if args.out_dir is not None else firedpy_parser.prompt_for_argument('out_dir')
    os.makedirs(out_dir, exist_ok=True)

    # Resolve full csv
    full_csv = str_to_bool(args.full_csv) if args.full_csv is not None else \
        firedpy_parser.prompt_for_argument('full_csv')

    n_cores = firedpy_parser.prompt_for_argument('n_cores') if args.n_cores is None else args.n_cores
    n_cores = os.cpu_count() - 1 if n_cores == 0 else n_cores

    # Resolve tile choice
    tile_choice = TileChoice(firedpy_parser.prompt_for_argument('tile_choice')) if args.tile_choice is None else (
        args.tile_choice)

    if tile_choice == TileChoice.A:
        tile_name = get_tile_name_from_directory(firedpy_parser, os.path.join(PROJECT_DIR, 'ref', 'continents'),
                                                 "Please enter the continent name:") if args.tile_name \
                                                                                        is None else args.tile_name
        shape_file = os.path.join("ref", "continents", tile_name + ".gpkg")
        print(f"Filtering for MODIS tiles that intersect \n  {shape_file}")
        tiles = shape_to_tiles(shape_file)
        eco_region_type, eco_region_level = (EcoRegionType.NA if tile_name == "north_america" else EcoRegionType.WORLD,
                                             3 if tile_name == "north_america" else None)

    elif tile_choice == TileChoice.B:
        tile_name = get_tile_name_from_directory(firedpy_parser,
                                                 os.path.join(PROJECT_DIR, 'ref', 'individual_countries'),
                                                 "Please enter the country name:") if args.tile_name \
                                                                                      is None else args.tile_name
        shape_file = os.path.join('ref', 'individual_countries', tile_name + ".gpkg")
        print(f"Filtering for MODIS tiles that intersect \n  {shape_file}")
        tiles = shape_to_tiles(shape_file)
        eco_region_type, eco_region_level = (EcoRegionType.NA, 3) if tile_name in [
            'united_States_of_america', 'canada', 'united_states_virgin_islands'] else (EcoRegionType.WORLD, None)

    elif tile_choice == TileChoice.C:
        tile_name = get_tile_name_from_directory(firedpy_parser, os.path.join(PROJECT_DIR, 'ref', 'us_states'),
                                                 "Please enter the state name:") if args.tile_name \
                                                                                    is None else args.tile_name
        shape_file = os.path.join('ref', 'us_states', tile_name + ".gpkg")
        print(f"Filtering for MODIS tiles that intersect \n  {shape_file}")
        tiles = shape_to_tiles(shape_file)
        eco_region_type, eco_region_level = EcoRegionType.NA, 3

    elif tile_choice == TileChoice.D:
        tiles = firedpy_parser.prompt_for_argument(
            arg_name='tile_name',
            prompt_override="Please enter tiles as a list of characters (no quotes no spaces)(e.g., h08v04 h09v04 " \
                            "...):").split(' ') if args.tile_name is None else args.tile_name
        tile_name = tiles[0]
        eco_region_type = firedpy_parser.prompt_for_argument('eco_region_type')
        eco_region_level = 3 if eco_region_type == EcoRegionType.NA else None
        shape_file = None
    else:
        raise ValueError('Invalid tile choice')

    daily = firedpy_parser.prompt_for_argument('daily') if args.daily is None else str_to_bool(args.daily)
    spatial_param = firedpy_parser.prompt_for_argument('spatial') if args.spatial_param is None else \
        args.spatial_param
    temporal_param = firedpy_parser.prompt_for_argument('temporal') if args.temporal_param is None else \
        args.temporal_param
    shape_type = ShapeType(
        firedpy_parser.prompt_for_argument('shape_type')) if args.shape_type is None else args.shape_type
    shapefile = shape_type != ShapeType.NONE

    # Resolve land cover type
    land_cover_type = (LandCoverType(firedpy_parser.prompt_for_argument('land_cover_type')) if args.land_cover_type is
                                                                                               None else
                       LandCoverType(args.land_cover_type))

    cleanup = str_to_bool(args.cleanup) if args.cleanup is not None else firedpy_parser.prompt_for_argument('cleanup')

    username = os.environ.get('FIREDPY_ED_USER', None)
    password = os.environ.get('FIREDPY_ED_PWD', None)
    if username is None or password is None:
        print("Please input your NASA Earthdata username and password in order to download the land cover data. If"
              " you do not have an Earthdata account, you can register at https://urs.earthdata.nasa.gov/. To "
              "avoid seeing this prompt again, you can set the FIREDPY_ED_USER and FIREDPY_ED_PWD environment "
              " variables.")
        username = firedpy_parser.prompt_for_argument('username')
        password = firedpy_parser.prompt_for_argument('password', sensitive=True)
        test_earthdata_credentials(username, password)

    start_year = firedpy_parser.prompt_for_argument('start_year') if args.start_year is None else args.start_year
    end_year = firedpy_parser.prompt_for_argument('end_year') if args.end_year is None else args.end_year
    start_year = None if start_year == 0 else start_year
    end_year = None if end_year == 0 else end_year

    if land_cover_type != LandCoverType.NONE:
        print('Retrieving landcover...')
        land_cover = LandCover(out_dir, n_cores=n_cores, username=username, password=password)
        land_cover.get_land_cover(tiles, land_cover_type)

    eco_region_data = EcoRegion(out_dir)
    eco_region_data.get_eco_region()

    burn_data = BurnData(out_dir, username, password, n_cores)
    burn_data.get_burns(tiles, start_year, end_year)

    # Create Model Builder object
    models = ModelBuilder(out_dir=out_dir, tiles=tiles, spatial_param=spatial_param, temporal_param=temporal_param,
                          n_cores=n_cores)

    event_perimeters = models.build_events()

    # TODO: This can be parallelized
    gdf = models.build_points(event_perimeters, shape_file_path=shape_file)
    gdf = models.add_fire_attributes(gdf)
    if land_cover_type != LandCoverType.NONE:
        gdf = models.add_land_cover_attributes(gdf, land_cover_type)
    gdf = models.process_geometry(gdf)
    gdf = models.add_eco_region_attributes(gdf, eco_region_type, eco_region_level)

    def generate_path(proj_dir, base_filename, shape_type: ShapeType):
        """Generate the appropriate file path."""
        file_ext = {
            ShapeType.SHP: ['.shp', None],
            ShapeType.GPKG: ['.gpkg', None],
            ShapeType.BOTH: ['.shp', '.gpkg']
        }.get(shape_type, ['.gpkg'])

        if isinstance(file_ext, list):
            return [os.path.join(proj_dir, "outputs", "shapefiles", base_filename + ext) if ext is not None else None
                    for ext in file_ext]

    # Build the polygons
    date_range = burn_data.get_date_range(start_year=start_year, end_year=end_year)
    base_file_name = f"fired_{tile_name}_{date_range[0][0]}_to_{date_range[-1][0]}"
    daily_base = f"{base_file_name}_daily"
    event_base = f"{base_file_name}_events"

    daily_shape_path, daily_gpkg_path = generate_path(out_dir, daily_base, shape_type)
    event_shape_path, event_gpkg_path = generate_path(out_dir, event_base, shape_type)
    csv_path = os.path.join(out_dir, base_file_name + '.csv')

    if daily:
        os.makedirs(os.path.dirname(daily_shape_path), exist_ok=True)
        if daily_gpkg_path is not None:
            os.makedirs(os.path.dirname(daily_gpkg_path), exist_ok=True)
        gdf = models.process_daily_data(gdf, csv_path.replace('.csv', '_daily.csv'), daily_shape_path, daily_gpkg_path)
    else:
        gdf = models.process_event_data(gdf)

    os.makedirs(os.path.dirname(event_shape_path), exist_ok=True)
    if event_gpkg_path is not None:
        os.makedirs(os.path.dirname(event_gpkg_path), exist_ok=True)
    models.save_event_data(gdf, csv_path, event_shape_path, event_gpkg_path, full_csv=full_csv)

    end = time.perf_counter()
    seconds = end - start
    minutes = seconds / 60
    print("Job completed in {} minutes".format(round(minutes, 2)))

    peak_mem = (peak_memory() - initial_memory) / (1024 * 1024 * 1024)
    print("Peak memory usage:", peak_mem, "GB")

    make_read_me(out_dir, tile_name, base_file_name, 1, date_range[0], date_range[-1],
                 daily, spatial_param, temporal_param, shapefile, shape_type, seconds, peak_mem, n_cores)

    if cleanup:
        cleanup_intermediate_files(out_dir)


if __name__ == "__main__":
    main()
