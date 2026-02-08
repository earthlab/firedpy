"""firedpy run methods."""
import os
import resource
import shutil
import urllib.request

from http.cookiejar import CookieJar
from logging import getLogger
from pathlib import Path

from firedpy.data_classes import ATTR_DESCS, BurnData, EcoRegion, LandCover
from firedpy.enums import ShapeType
from firedpy.model_classes import ModelBuilder
from firedpy.utilities.create_readme import make_read_me
from firedpy.utilities.logging import init_logger
from firedpy.utilities.spatial import country_to_tiles

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


def get_peak_memory():
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


def fired(
    project_directory,
    project_name=None,  # Changing this soon to like project name
    country="Iceland",
    tiles=None,
    shape_file=None,
    start_year=2000,
    end_year=2025,
    spatial_param=8,
    temporal_param=3,
    daily=True,
    shape_type="gpkg",
    eco_region_level=1,
    eco_region_type=None,
    land_cover_type=None,
    full_csv=True,
    n_cores=0,
    cleanup=False
):
    """Run all steps of the firedpy modeling pipeline.

    TODO: Consolidate some of these code chunks to make it a bit easier to
        read.

    Parameters
    ----------
    project_directory : str
        Project output directory path. Required.
    project_name : str | NoneType
        A name used to identify the output files of this project. Defaults to
        None, which will use the name of the parent run directory.
    country : str
        The name of a country to use as a study area. Defaults to 'Iceland'.
    tiles : str | list
        A string representing a single MODIS tile (e.g., 'h08v04'), a string
        representing multiple tiles separated by spaces (e.g., 'h08v04 h09v04')
        or a list representing multiple tiles (e.g., ['h08v04', 'h09v04']).
        If None, a `country` or `shape_file` parameter is required.
    shape_file : str
        Path to a shapefile to use for the fire study area. Defaults to None.
        If not provided, a `tiles` or `country` parameter is required.
    start_year : int
        The first year of fire events. Defaults to 2000.
    end_year : int
        The last year of fire events. Defaults to 2025.
    spatial_param : int
        The number of cells (~463 m resolution) to search for neighboring burn
        detections. Defaults to 5 cells in all directions.
    temporal_param : int
        The number of days to search for neighboring burn detections.
    daily : boolean
        Create the daily polygons or just the event-level perimeter for your
        analysis area. If this flag is set, the daily and event polygons will
        be created, otherwise only the event level.
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
            4: Annual BIOME-Biogeochemical Cycles (BGC)
            5: Annual Plant Functional Types (PFT)

        If you do not have an account register at
        https://urs.earthdata.nasa.gov/home.

        Defaults to None.
    full_csv : bool
        Export all attributes to CSV. Defaults to only x and y coordinates,
        event date, and event id will be exported to a CSV. Defaults to
        True.
    n_cores : int
        Number of cores to use for parallel processing. Defaults to 0
        or all available cores.
    cleanup : bool
        Cleanup. If set then the burn area and landcover files will be removed
        after each run to save disk space in between multiple runs. Defaults
        to False.

    Returns
    -------
    geopandas.geodataframe.GeoDataFrame : A geodataframe of fire event
        perimeters and attributes.
    """
    # Setup logging for this output directory
    out_dir = Path(project_directory).expanduser().absolute()
    init_logger(out_dir=out_dir)
    logger.info(
        f"Running firedpy for years {start_year} to {end_year} on MODIS "
        f"tiles: {tiles}."
    )

    # Format study area parameters
    if country:
        tiles = country_to_tiles(country)
    else:
        if isinstance(tiles, str):
            tiles = tiles.replace("'", "").split()

    # Get the burn data
    logger.info("Collecting MODIS burn data.")
    burn_data = BurnData(out_dir=out_dir, n_cores=n_cores)
    burn_data.get_burns(
        tiles=tiles,
        start_year=start_year,
        end_year=end_year
    )

    # Create a model builder object
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

    # Add initial fire characteristic
    gdf = models.add_fire_attributes(gdf=gdf)

    # Add land cover characteristic if requested
    if land_cover_type:
        lc_desc = ATTR_DESCS["landcover"][str(land_cover_type)]
        logger.info(f"Adding landcover attributes: {lc_desc}.")
        land_cover = LandCover(
            out_dir=out_dir,
            n_cores=n_cores
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
    if eco_region_type:
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

    # Set output paths
    out_dir = os.path.expanduser(out_dir)
    date_range = burn_data.get_date_range(
        start_year=start_year,
        end_year=end_year
    )
    date1 = date_range[0][0]
    date2 = date_range[-1][0]
    if not project_name:
        project_name = Path(project_directory).name
    base_file_name = f"fired_{project_name}_{date1}_to_{date2}"
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

    model_outputs_dir = Path(out_dir).joinpath("outputs")
    shape_dir = model_outputs_dir.joinpath("shapefiles")
    table_dir = model_outputs_dir.joinpath("tables")
    csv_path = table_dir.joinpath(f"{base_file_name}.csv")
    shape_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(exist_ok=True)

    # Process event data
    if daily:
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

    # make_read_me(
    #     gdf=gdf,
    #     project_directory=project_directory,
    #     tiles=tiles,
    #     spatial_param=spatial_param,
    #     temporal_param=temporal_param,
    #     shapefile=shape_file,
    #     runtime=runtime,
    #     n_cores=n_cores,
    #     peak_memory=peak_memory,
    # )

    # # Remove intermediate files if requested
    # if cleanup:
    #     cleanup_intermediate_files(project_directory)

    return gdf


if __name__ == "__main__":
    project_directory = '/home/travis/scratch/firedpy/test2'
    interactive = True
    tiles = None
    country = 'United States of America'
    shape_file = None
    start_year = 2000
    end_year = 2025
    daily = False
    spatial_param = 8  # pixels (nominally ~3,704 m but varies by location)
    temporal_param = 3  # days
    shape_type = 'gpkg'  # GeoPackage
    eco_region_level = 1  # Level I - Least Detailed
    eco_region_type = 'na'  # North American Ecoregions (Omernick, 1987)
    land_cover_type = 1  # International Geosphere-Biosphere Programme (IGBP) scheme
    full_csv = False
    n_cores = 32
    cleanup = False
