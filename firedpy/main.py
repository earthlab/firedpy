import argparse
import os
import resource
import shutil
import time
import urllib.request
import warnings

from http.cookiejar import CookieJar
from pathlib import Path

from firedpy import DATA_DIR
from firedpy.data_classes import BurnData, EcoRegion, LandCover
from firedpy.enums import (
    EcoRegionType,
    LandCoverType,
    ShapeType,
    TileChoice
)
from firedpy.help import HELP_TEXT
from firedpy.model_classes import ModelBuilder
from firedpy.utilities.argument_parser import FiredpyArgumentParser
from firedpy.utilities.create_readme import make_read_me
from firedpy.spatial import shape_to_tiles

warnings.filterwarnings("ignore", category=FutureWarning)

EARTHDATA_TEST_URL = (
    "https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.061/2019.01.01/"
    "BROWSE.MCD12Q1.A2019001.h10v09.061.2022169160720.1.jpg"
)


def generate_path(proj_dir, base_filename, shape_type: ShapeType):
    """Generate the appropriate file path."""
    proj_dir = Path(proj_dir)
    file_exts = {
        ShapeType.SHP: [".shp", None],
        ShapeType.GPKG: [".gpkg", None],
        ShapeType.BOTH: [".shp", ".gpkg"]
    }
    file_ext = file_exts.get(shape_type, [".gpkg"])

    paths = []
    for ext in file_ext:
        if ext:
            fname = f"{base_filename}{ext}"
            path = proj_dir.joinpath("outputs", "shapefiles", fname)
        else:
            path = None
        paths.append(path)

    return paths


def get_tile_name_from_directory(parser: FiredpyArgumentParser, directory: str,
                                 prompt_message: str):
    """Function to prompt the user for a tile name from a given directory."""
    available_files = [
        file.replace(".gpkg", "") for file in os.listdir(directory) if
        file.endswith(".gpkg")
    ]
    return parser.prompt_for_argument(
        arg_name="tile_name",
        prompt_override=prompt_message,
        accepted_value_override=available_files
    )


def str_to_bool(s: str):
    accepted_values = ["true", "y", "yes", "false", "n", "no"]
    if s.lower() not in accepted_values:
        raise ValueError(f"{s} must be in {accepted_values}")
    return s.lower() in ["true", "y", "yes"]


def test_earthdata_credentials(username: str, password: str) -> None:
    # Earthdata Login
    password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    password_manager.add_password(None, "https://urs.earthdata.nasa.gov",
                                  username, password)

    # Create a cookie jar for storing cookies. This is used to store and return
    # the session cookie given to use by the data server (otherwise it will
    # just keep sending us back to Earthdata Login to authenticate).  Ideally,
    # we should use a file based cookie jar to preserve cookies between runs.
    # This will make it much more efficient.
    cookie_jar = CookieJar()

    # Install all the handlers
    opener = urllib.request.build_opener(
        urllib.request.HTTPBasicAuthHandler(password_manager),
        # urllib.request.HTTPHandler(debuglevel=1),  # Uncomment to see details
        # urllib.request.HTTPSHandler(debuglevel=1),  # of requests/responses
        urllib.request.HTTPCookieProcessor(cookie_jar))
    urllib.request.install_opener(opener)

    # Send a test URL
    request = urllib.request.Request(EARTHDATA_TEST_URL)
    urllib.request.urlopen(request)


def peak_memory():
    """Get maximum resident usage size of current process."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


def cleanup_intermediate_files(out_dir):
    shutil.rmtree(os.path.join(out_dir, "rasters", "burn_area"))
    shutil.rmtree(os.path.join(out_dir, "rasters", "land_cover"))


def main():
    # Start the timer (seconds)
    start = time.perf_counter()

    # Get the maximum resource use for this process
    initial_memory = peak_memory()

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--output_directory",
        dest="out_dir",
        default=DATA_DIR.joinpath("output"),
        help=HELP_TEXT["data"]
    )
    parser.add_argument(
        "-f",
        "--file_name",
        dest="file_name",
        default="fired",
        type=str,
        help=HELP_TEXT["file"]
    )
    parser.add_argument(
        "-et",
        "--eco_region_type",  # Distinguish this one from the eco region level
        type=EcoRegionType,
        help=HELP_TEXT["eco"]
    )
    parser.add_argument(
        "-el",
        "--eco_region_level",
        type=int,
        help=HELP_TEXT["eco"]
    )
    parser.add_argument(
        "-lc",
        "--land_cover_type",
        help=HELP_TEXT["lc"],
        type=int
    )
    parser.add_argument(
        "-st",
        "--shape_type",
        help=HELP_TEXT["shp"],
        type=ShapeType
    )
    parser.add_argument(
        "-s",
        "--spatial",
        dest="spatial_param",
        type=int,
        help=HELP_TEXT["sp"]
    )
    parser.add_argument(
        "-t",
        "--temporal",
        dest="temporal_param",
        type=int,
        help=HELP_TEXT["tmp"]
    )
    parser.add_argument(
        "-tc",
        "--tile_choice",
        type=TileChoice,
        help=""
    )
    parser.add_argument(
        "-tn",
        "--tile_name",
        type=str,
        help=HELP_TEXT["tile_name"]
    )
    parser.add_argument(
        "-d",
        "--daily",
        type=str,
        help=HELP_TEXT["daily"]
    )
    parser.add_argument(
        "-sy",
        "--start_year",
        type=int,
        help=HELP_TEXT["start_year"]
    )
    parser.add_argument(
        "-ey",
        "--end_year",
        type=int,
        help=HELP_TEXT["end_year"]
    )
    parser.add_argument(
        "-fc",
        "--full_csv",
        type=str,
        help=HELP_TEXT["full_csv"]
    )
    parser.add_argument(
        "-nc",
        "--n_cores",
        type=int,
        help=HELP_TEXT["n_cores"]
    )
    parser.add_argument(
        "-cu",
        "--cleanup",
        type=str,
        help=HELP_TEXT["cleanup"]
    )
    args = parser.parse_args()

    param_fpath = DATA_DIR.joinpath("params.txt")
    firedpy_parser = FiredpyArgumentParser(param_fpath)

    # Resolve the output directory
    if args.out_dir:
        out_dir = args.out_dir
    else:
        out_dir = firedpy_parser.prompt_for_argument("out_dir")
    os.makedirs(out_dir, exist_ok=True)

    # Resolve full csv
    if args.full_csv:
        full_csv = str_to_bool(args.full_csv)
    else:
        firedpy_parser.prompt_for_argument("full_csv")

    # Resolve n_cores
    if args.n_cores:
        n_cores = firedpy_parser.prompt_for_argument("n_cores")
    else:
        n_cores = args.n_cores
    if n_cores == 0:
        n_cores = os.cpu_count() - 1

    # Resolve tile choice
    if args.tile_choice:
        tile_choice = args.tile_choice
    else:
        user_input = firedpy_parser.prompt_for_argument("tile_choice")
        tile_choice = TileChoice(user_input)

    if tile_choice == TileChoice.A:
        if args.tile_name:
            tile_name = args.tile_name
        else:
            tile_name = get_tile_name_from_directory(
                parser=firedpy_parser,
                directory=DATA_DIR.joinpath("continents"),
                prompt_message="Please enter the continent name: "
            )

        shape_file = DATA_DIR.joinpath("continents", f"{tile_name}.gpkg")
        print(f"Filtering for MODIS tiles that intersect \n  {shape_file}")
        tiles = shape_to_tiles(shape_file)

        if tile_name == "north_america":
            eco_region_type = EcoRegionType.NA
            eco_region_level = 3
        else:
            eco_region_type = EcoRegionType.WORLD
            eco_region_level = None

    elif tile_choice == TileChoice.B:
        if args.tile_name:
            tile_name = args.tile_name
        else:
            tile_name = get_tile_name_from_directory(
                parser=firedpy_parser,
                directory=DATA_DIR.joinpath("individual_countries"),
                prompt_message="Please enter the country name:"
            )

        shape_file = DATA_DIR.joinpath("individual_countries",
                                       f"{tile_name}.gpkg")
        print(f"Filtering for MODIS tiles that intersect \n  {shape_file}")
        tiles = shape_to_tiles(shape_file)

        na_tiles = ["united_states_of_america", "canada",
                    "united_states_virgin_islands"]
        if tile_name in na_tiles:
            eco_region_type = EcoRegionType.NA
            eco_region_level = 3
        else:
            eco_region_type = EcoRegionType.WORLD
            eco_region_level = None

    elif tile_choice == TileChoice.C:
        if args.tile_name:
            tile_name = args.tile_name
        else:
            tile_name = get_tile_name_from_directory(
                parser=firedpy_parser,
                directory=DATA_DIR.joinpath("us_states"),
                prompt_message="Please enter the state name:"
            )

        shape_file = DATA_DIR.joinpath("us_states", f"{tile_name}.gpkg")
        print(f"Filtering for MODIS tiles that intersect \n  {shape_file}")
        tiles = shape_to_tiles(shape_file)
        eco_region_type = EcoRegionType.NA
        eco_region_level = 3

    elif tile_choice == TileChoice.D:
        shape_file = input("Please enter the GPKG or shapefile path:")
        print(f"Filtering for MODIS tiles that intersect \n  {shape_file}")
        tiles = shape_to_tiles(shape_file)
        tile_name = os.path.basename(shape_file).rstrip(".gpkg").rstrip(".shp")
        eco_region_type = firedpy_parser.prompt_for_argument("eco_region_type")
        if eco_region_type == EcoRegionType.NA:
            eco_region_level = 3
        else:
            eco_region_level = None

    elif tile_choice == TileChoice.E:
        tiles = firedpy_parser.prompt_for_argument(
            arg_name="tile_name",
            prompt_override=(
                "Please enter tiles as a list of characters "
                "(no quotes no spaces)(e.g., h08v04 h09v04 ...):"
            )
        )
        if args.tile_name:
            tiles = args.tile_name
        else:
            tiles = tiles.split(" ")
        tile_name = tiles[0]
        eco_region_type = firedpy_parser.prompt_for_argument("eco_region_type")
        if eco_region_type == EcoRegionType.NA:
            eco_region_level = 3
        else:
            eco_region_level = None
        shape_file = None
    else:
        raise ValueError("Invalid tile choice.")

    # Resolve daily option
    if args.daily:
        daily = str_to_bool(args.daily)
    else:
        daily = firedpy_parser.prompt_for_argument("daily")

    # Resolve spatial parameter
    if args.spatial_param:
        spatial_param = args.spatial_param
    else:
        spatial_param = firedpy_parser.prompt_for_argument("spatial")

    # Resolve temporal parameter
    if args.temporal_param:
        temporal_param = args.temporal_param
    else:
        temporal_param = firedpy_parser.prompt_for_argument("temporal")

    # Resolve shapefile format option
    if args.shape_type:
        shape_type = args.shape_type
    else:
        user_input = firedpy_parser.prompt_for_argument("shape_type")
        shape_type = ShapeType(user_input)
    shapefile = shape_type != ShapeType.NONE

    # Resolve land cover type
    if args.land_cover_type:
        land_cover_type = LandCoverType(args.land_cover_type)
    else:
        user_input = firedpy_parser.prompt_for_argument("land_cover_type")
        land_cover_type = LandCoverType(user_input)

    # Resolve clean up option
    if args.cleanup:
        cleanup = str_to_bool(args.cleanup)
    else:
        cleanup = firedpy_parser.prompt_for_argument("cleanup")

    # Resolve username and password
    username = os.environ.get("FIREDPY_ED_USER", None)
    password = os.environ.get("FIREDPY_ED_PWD", None)
    if username is None or password is None:
        print(HELP_TEXT["earthdata"])
        username = firedpy_parser.prompt_for_argument("username")
        password = firedpy_parser.prompt_for_argument(
            "password",
            sensitive=True
        )
        print("EarthAccess will handle authentication automatically.")

    # Resolve start and end years
    if args.start_year and args.start_year > 0:
        start_year = args.start_year
    elif args.start_year == 0:
        start_year = None
    else:
        start_year = firedpy_parser.prompt_for_argument("start_year")

    if args.end_year and args.end_year > 0:
        end_year = args.end_year
    elif args.end_year == 0:
        end_year = None
    else:
        end_year = firedpy_parser.prompt_for_argument("end_year")

    # Resolve land cover
    if land_cover_type != LandCoverType.NONE:
        print("Retrieving landcover...")
        land_cover = LandCover(out_dir, n_cores=n_cores, username=username,
                               password=password)
        land_cover.get_land_cover(tiles, land_cover_type)

    # Get the ecoregion data
    eco_region_data = EcoRegion(out_dir)
    eco_region_data.get_eco_region()

    # Get the burn data
    burn_data = BurnData(out_dir, username, password, n_cores)
    burn_data.get_burns(tiles, start_year, end_year)

    # Create Model Builder object
    models = ModelBuilder(
        out_dir=out_dir,
        tiles=tiles,
        spatial_param=spatial_param,
        temporal_param=temporal_param,
        n_cores=n_cores
    )

    # Build event perimeters
    event_perimeters = models.build_events()

    # TODO: This can be parallelized
    gdf = models.build_points(event_perimeters, shape_file_path=shape_file)
    gdf = models.add_fire_attributes(gdf)
    if land_cover_type != LandCoverType.NONE:
        gdf = models.add_land_cover_attributes(gdf, tiles, land_cover_type)
    gdf = models.process_geometry(gdf)
    gdf = models.add_eco_region_attributes(gdf, eco_region_type,
                                           eco_region_level)

    # Calculate fire spread speed and maximum travel vectors
    gdf = models.add_kg_attributes(gdf)

    # Build the polygons
    date_range = burn_data.get_date_range(
        start_year=start_year,
        end_year=end_year
    )
    date1 = date_range[0][0]
    date2 = date_range[-1][0]
    base_file_name = f"fired_{tile_name}_{date1}_to_{date2}"
    daily_base = f"{base_file_name}_daily"
    event_base = f"{base_file_name}_events"
    daily_shape_path, daily_gpkg_path = generate_path(
        proj_dir=out_dir,
        base_filename=daily_base,
        shape_type=shape_type
    )
    event_shape_path, event_gpkg_path = generate_path(
        proj_dir=out_dir,
        base_filename=event_base,
        shape_type=shape_type
    )
    csv_path = os.path.join(out_dir, base_file_name + ".csv")

    if daily:
        os.makedirs(os.path.dirname(daily_shape_path), exist_ok=True)
        csv_path = str(csv_path).replace(".csv", "_daily.csv")
        if daily_gpkg_path:
            os.makedirs(os.path.dirname(daily_gpkg_path), exist_ok=True)
        gdf = models.process_daily_data(gdf, csv_path, daily_shape_path,
                                        daily_gpkg_path)
    else:
        gdf = models.process_event_data(gdf)

    os.makedirs(os.path.dirname(event_shape_path), exist_ok=True)
    if event_gpkg_path is not None:
        os.makedirs(os.path.dirname(event_gpkg_path), exist_ok=True)
    models.save_event_data(gdf, csv_path, event_shape_path, event_gpkg_path,
                           full_csv=full_csv)

    end = time.perf_counter()
    seconds = end - start
    minutes = seconds / 60
    print("Job completed in {} minutes".format(round(minutes, 2)))

    peak_mem = (peak_memory() - initial_memory) / (1024 * 1024 * 1024)
    print("Peak memory usage:", peak_mem, "GB")

    make_read_me(out_dir, tile_name, base_file_name, 1, date_range[0],
                 date_range[-1], daily, spatial_param, temporal_param,
                 shapefile, shape_type, seconds, peak_mem, n_cores)

    if cleanup:
        cleanup_intermediate_files(out_dir)


if __name__ == "__main__":
    main()
