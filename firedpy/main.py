import os
import resource
import shutil
import time
import urllib.request

from http.cookiejar import CookieJar
from logging import getLogger
from pathlib import Path

from firedpy.data_classes import ATTR_DESCS, BurnData, EcoRegion, LandCover
from firedpy.enums import LandCoverType, ShapeType
from firedpy.model_classes import ModelBuilder
# from firedpy.utilities.argument_parsing import FiredpyArgumentParser
from firedpy.utilities.create_readme import make_read_me
from firedpy.utilities.logging import init_logger

logger = getLogger(__name__)


EARTHDATA_TEST_URL = (
    "https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.061/2019.01.01/"
    "BROWSE.MCD12Q1.A2019001.h10v09.061.2022169160720.1.jpg"
)


def cleanup_intermediate_files(out_dir):
    """Remove temporary `burn_area` and `land_cover` raster files.

    Parameters
    ----------
    out_dir : str | pathlib.PosixPath
        The output directory containing 'rasters/burn_area' and
        'rasters/land_cover'.
    """
    shutil.rmtree(os.path.join(out_dir, "rasters", "burn_area"))
    shutil.rmtree(os.path.join(out_dir, "rasters", "land_cover"))


def generate_path(proj_dir, base_filename, shape_type):
    """Return the full paths to a target files.

    Parameters
    ----------
    proj_dir : str
        A project directory selected containing target shapefiles. This
        corresponds with the '-o' or '--output_directory` option in the
        firedpy CLI.
    base_filename : str
        The file name of the target file.
    shape_type : firedpy.enums.ShapeType
        One of ShapeType.SHP, ShapeType.GPKG, or ShapeType.BOTH.

    Returns
    -------
    list[str] : A list of full filepaths corresponding with the target file
        and file formats.
    """
    # What's the benefit of enums here?
    if not isinstance(shape_type, ShapeType):
        shape_type = ShapeType(shape_type)

    # Get the requested combination of shapefile types
    file_extensions = {
        ShapeType.SHP: [".shp", None],
        ShapeType.GPKG: [None, ".gpkg"],
        ShapeType.BOTH: [".shp", ".gpkg"]
    }
    file_ext = file_extensions.get(shape_type, [".gpkg"])

    # Create the paths
    paths = []
    proj_dir = os.path.expanduser(proj_dir)
    for ext in file_ext:
        if ext:
            fname = f"{base_filename}{ext}"
            path = os.path.join(proj_dir, "outputs", "shapefiles", fname)
        else:
            path = None
        paths.append(path)

    return paths


def peak_memory():
    """Get maximum resident usage size of current process.

    NOTE: This probably isn't doing what we want.
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


def test_earthdata_credentials(username, password):
    """Test access to `ers.earthdata.nasa.gov`

    Parameters
    ----------
    username : str
        Username for Earthdata account.
    password : str
        Password for Earthdata account.
    """
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


def run_all(
    out_dir,
    tiles,
    tile_name=None,
    start_year=2000,
    end_year=2025,
    daily=True,
    spatial_param=8,
    temporal_param=3,
    shape_file=None,
    shape_type="gpkg",
    eco_region_level=None,
    eco_region_type="na",
    land_cover_type=1,
    n_cores=0,
    full_csv=True,
    username=None,
    password=None
):
    """Run all steps of the firedpy modeling pipeline.

    Parameters
    ----------
    out_dir : str
        Project output directory path. Required.
    tiles : list
        List of MODIS tiles (e.g., ['h08v04', 'h09v04']). Required.
    tile_name : str | NoneType
        The name of the MODIS tile being run? Shouldn't there be multiple?
        How is this different from above, defaulting to None for now.
    start_year : int
        The first year of fire events. Defaults to 2000.
    end_year : int
        The last year of fire events. Defaults to 2025.
    daily : boolean
        Create the daily polygons or just the event-level perimeter for your
        analysis area. If this flag is set, the daily and event polygons will
        be created, otherwise only the event level.
    spatial_param : int
        The number of cells (~463 m resolution) to search for neighboring burn
        detections. Defaults to 5 cells in all directions.
    temporal_param : int
        The number of days to search for neighboring burn detections.
    shape_file : str
        Path to a shapefile to use for the fire study area.
    shape_type : str
        Build shapefiles from the event data frame. Specify either "shp",
        "gpkg", or both. Shapefiles of both daily progression and overall
        event perimeters will be written to the 'outputs/shapefiles' folder of
        the chosen project directory. These will be saved in the specified
        geopackage format (.gpkg), ERSI Shapefile format (.shp), or save them
        in both formats using the file basename of the fire event data frame
        (e.g. 'modis_events_daily.gpkg' and 'modis_events.gpkg').
    eco_region_level : int
        The desired Ecoregions level from the North American Commission for
        Environmental Cooperation (CEC). Levels 1 to 3 are available, with
        level 1 representing the broadest scale and level III representing the
        most detailed. Defaults to 1.
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
    land_cover_type : int
        Include land cover as an attribute, provide a number corresponding
        with a MODIS/Terra+Aqua Land Cover (MCD12Q1) category followed with
        username:password of your NASA's Earthdata service account.
        Available land cover categories:

            1: IGBP global vegetation classification scheme
            2: University of Maryland (UMD) scheme
            3: MODIS-derived LAI/fPAR scheme

        If you do not have an account register at
        https://urs.earthdata.nasa.gov/home. Defaults to 1.
    n_cores : int
        Number of cores to use for parallel processing. Defaults to 0
        or all available cores.
    full_csv : bool
        Export all attributes to CSV. Defaults to only x and y coordinates,
        event date, and event id will be exported to a CSV. Defaults to
        True.
    username : str
        Username for a NASA Earthdata Account. Defaults to None, will
        prompt user.
    password : str
        Password for a NASA Earthdata Account. Defaults to None, will
        prompt user.
    """
    # Setup logging for this output directory
    out_dir = Path(out_dir).expanduser()
    init_logger(out_dir=out_dir)
    logger.info(
        f"Running firedpy for years {start_year} to {end_year} on MODIS "
        f"tiles: {tiles}."
    )

    # Get the burn data
    logger.info("Collecting MODIS burn data.")
    burn_data = BurnData(
        out_dir=out_dir,
        username=username,
        password=password,
        n_cores=n_cores
    )
    burn_data.get_burns(
        tiles=tiles,
        start_year=start_year,
        end_year=end_year
    )

    # Create Model Builder object
    models = ModelBuilder(
        out_dir=out_dir,
        tiles=tiles,
        spatial_param=spatial_param,
        temporal_param=temporal_param,
        n_cores=n_cores
    )

    # Build event perimeters
    logger.info("Building event perimeter geometries.")
    event_perimeters = models.build_events()

    # Build the event geodataframe
    gdf = models.build_points(
        event_perimeters=event_perimeters,
        shape_file=shape_file
    )

    # Add secondary attributes
    gdf = models.add_fire_attributes(gdf=gdf)
    if land_cover_type != LandCoverType.NONE:
        lc_desc = ATTR_DESCS["landcover"][land_cover_type]
        logger.info(f"Adding landcover attributes: {lc_desc}.")
        land_cover = LandCover(
            out_dir=out_dir,
            n_cores=n_cores,
            username=username,
            password=password
        )
        land_cover.get_land_cover(
            tiles=tiles,
            land_cover_type=land_cover_type
        )
        gdf = models.add_land_cover_attributes(
            gdf=gdf,
            tiles=tiles,
            land_cover_type=land_cover_type
        )

    # What is this doing? (buffering with an envelope)
    gdf = models.process_geometry(gdf)

    # Add ecoregion attributes
    eco_desc = ATTR_DESCS["ecoregions"][eco_region_type]
    msg = f"Adding ecoregion data: {eco_desc}, Level {eco_region_level}."
    logger.info(msg)
    eco_region_data = EcoRegion(out_dir=out_dir)
    eco_region_data.get_eco_region()
    gdf = models.add_eco_region_attributes(
        gdf=gdf,
        eco_region_type=eco_region_type,
        eco_region_level=eco_region_level
    )

    # Calculate fire spread speed and maximum travel vectors
    gdf = models.add_kg_attributes(gdf)

    # Set paths
    out_dir = os.path.expanduser(out_dir)
    date_range = burn_data.get_date_range(
        start_year=start_year,
        end_year=end_year
    )
    date1 = date_range[0][0]
    date2 = date_range[-1][0]
    if tile_name:
        base_file_name = f"fired_{tile_name}_{date1}_to_{date2}"
    else:
        base_file_name = f"fired_{date1}_to_{date2}"
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
    shape_dir = os.path.join(out_dir, "outputs", "shapefiles")
    table_dir = os.path.join(out_dir, "tables")
    csv_path = os.path.join(table_dir, base_file_name + ".csv")

    # Process event data
    if daily:
        os.makedirs(shape_dir, exist_ok=True)
        os.makedirs(table_dir, exist_ok=True)
        output_csv_path = str(csv_path).replace(".csv", "_daily.csv")
        gdf = models.process_daily_data(
            gdf=gdf,
            output_csv_path=output_csv_path,
            daily_shape_path=daily_shape_path,
            daily_gpkg_path=daily_gpkg_path
        )
    else:
        gdf = models.process_event_data(gdf)

    # Save the events
    models.save_event_data(
        gdf=gdf,
        output_csv_path=csv_path,
        event_shape_path=event_shape_path,
        event_gpkg_path=event_gpkg_path,
        full_csv=full_csv
    )

    return base_file_name


def main():
    # Start the timer (seconds, not as helpful for prompted inputs)
    start = time.perf_counter()

    # Get the maximum resource use for this process (fix this, it's broken)
    initial_memory = peak_memory()

    # Parse user arguments
    resolver = FiredpyArgumentParser()
    args = resolver.resolve()

    # Break out arguments for easier editing/readability
    out_dir = args.out_dir
    tile_name = args.tile_name
    tiles = args.tiles
    daily = args.daily
    start_year = args.date_range[0]
    end_year = args.date_range[1]
    spatial_param = args.spatial_param
    temporal_param = args.temporal_param
    shapefile = args.shapefile
    shape_type = args.shape_typem
    n_cores = args.n_cores
    full_csv = args.full_csv
    username = args.username
    password = args.password

    # Run the model
    file_base = run_all(
        out_dir=out_dir,
        tiles=tiles,
        tile_name=tile_name,
        start_year=start_year,
        end__year=end_year,
        daily=daily,
        spatial_param=spatial_param,
        temporal_param=temporal_param,
        shapefile=shapefile,
        shape_type=shape_type,
        n_cores=n_cores,
        full_csv=full_csv,
        username=username,
        password=password
    )

    # Done.
    end = time.perf_counter()
    seconds = end - start
    minutes = seconds / 60
    print(f"Job completed in {minutes:.2f} minutes")
    peak_mem = (peak_memory() - initial_memory) / 1_024 ** 3
    print(f"Peak memory usage: {peak_mem:.2f} GB")
    make_read_me(
        out_dir=out_dir,
        tile_name=tile_name,
        file_base=file_base,
        input=1,
        first_date=start_year,
        last_date=end_year,
        daily=daily,
        spatial_param=spatial_param,
        temporal=temporal_param,
        shapefile=shapefile,
        shape_type=shape_type,
        job_time=seconds,
        job_memory=peak_mem,
        n_cores=n_cores
    )

    # Remove intermediate files if requested
    if args.cleanup:
        cleanup_intermediate_files(out_dir)


if __name__ == "__main__":
    out_dir = "~/scratch/firedpy/test"
    tiles = ["h08v04", "h09v04"]
    tile_name = None
    start_year = 2020
    end_year = 2021
    daily = True
    spatial_param = 8
    temporal_param = 3
    shape_file = None
    shape_type = "gpkg"
    eco_region_level = 1
    eco_region_type = "na"
    land_cover_type = 1
    n_cores = 1
    full_csv = True
    username = None
    password = None

    key_fpath = os.path.expanduser("~/.keys/earthdata")
    if os.path.exists(key_fpath):
        with open(key_fpath) as file:
            lines = file.readlines()
        username = lines[0].strip()
        password = lines[1].strip()

    run_all(
        out_dir=out_dir,
        tiles=tiles,
        tile_name=tile_name,
        start_year=start_year,
        end_year=end_year,
        daily=daily,
        spatial_param=spatial_param,
        temporal_param=temporal_param,
        shape_file=shape_file,
        shape_type=shape_type,
        eco_region_level=eco_region_level,
        eco_region_type=eco_region_type,
        land_cover_type=land_cover_type,
        n_cores=n_cores,
        full_csv=full_csv,
        username=username,
        password=password
    )
