"""firedpy run methods."""
import os
import shutil
import time
import tracemalloc

from logging import getLogger
from pathlib import Path

from firedpy.data_classes import ATTR_DESCS, BurnData, EcoRegion, LandCover
from firedpy.model_classes import ModelBuilder
from firedpy.utilities.create_readme import make_read_me
from firedpy.utilities.logger import init_logger

logger = getLogger(__name__)


def cleanup_intermediate_files(project_directory):
    """Remove temporary  `burn_area`,  `land_cover`, and `eco_region` and raster files.

    Parameters
    ----------
    project_directory : str | pathlib.PosixPath
        The output directory containing 'rasters/burn_area' and
        'rasters/land_cover'.
    """
    shutil.rmtree(os.path.join(project_directory, "rasters", "burn_area"))
    shutil.rmtree(os.path.join(project_directory, "rasters", "land_cover"))

def fired(
    project_directory,
    project_name=None,
    country=None,
    tiles=None,
    shape_file=None,
    start_year=2000,
    end_year=2026,
    spatial_param=8,
    temporal_param=3,
    daily=True,
    shape_type="gpkg",
    eco_region_type=None,
    eco_region_level=3,
    land_cover_type=None,
    csv_type="none",
    n_cores=0,
    cleanup=False
):
    """Run all steps of the firedpy modeling pipeline.

    Parameters
    ----------
    project_directory : str
        Project output directory path. Required.
    project_name : str | NoneType
        A name used to identify the output files of this project. Defaults to
        None, which will use the name of the parent run directory.
    country : str
        The name of a country to use as a study area. If not provided, either
        'tiles' or 'shape_file' parameter must be provided.
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
        The last year of fire events. Defaults to 2026.
    spatial_param : int
        The number of cells (~463 m resolution) to search for neighboring burn
        detections.
    temporal_param : int
        The number of days to search for neighboring burn detections.
    daily : boolean
        Create the daily polygons or just the event-level perimeter for your
        analysis area. If this flag is set, the daily and event polygons will
        be created, otherwise only the event level.
    shape_type : str
        Build shapefiles from the event data frame. Specify either "shp",
        "gpkg", or both. Output files are written directly to the 'outputs/'
        folder of the chosen project directory in the specified geopackage
        format (.gpkg), ESRI Shapefile format (.shp), or both (e.g.
        'fired_events_daily.gpkg' and 'fired_events.gpkg').
    eco_region_level : int
        The desired Ecoregions level from the North American Commission for
        Environmental Cooperation (CEC). Levels 1 to 3 are available, with
        level 1 representing the broadest scale and level III representing the
        most detailed. Defaults to 0 (no eco region).
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
    csv_type : str
        Controls CSV output (case-insensitive). Options: 'full' (all
        attributes), 'events' (x, y, id, ig_date, last_date only), or 'none'
        (no CSV written). Defaults to 'none'.
    n_cores : int
        Number of cores to use for parallel processing. A value of 0 or None
        will use all available cores. Defaults to 0.
    cleanup : bool
        Cleanup. If set then the burn area and landcover files will be removed
        after each run to save disk space in between multiple runs. Defaults
        to False.

    Returns
    -------
    geopandas.geodataframe.GeoDataFrame : A geodataframe of fire event
        perimeters and attributes.
    """
    # Start the timer and memory tracer
    start = time.perf_counter()
    tracemalloc.start()

    # Setup logging for this output directory
    project_directory = Path(project_directory).expanduser().absolute()
    run_name = f"fired_{project_name}" if project_name else None
    init_logger(project_directory=project_directory, run_name=run_name)
    logger.info(
        f"Running firedpy for years {start_year} to {end_year} on MODIS "
        f"tiles: {tiles}."
    )

    # Make sure they have a target location
    if not country and not tiles and not shape_file:
        msg = ("Must provide one of `country`, `tiles`, or `shape_file` "
               "parameters")
        logger.error(msg)
        raise KeyError(msg)

    # Format study area parameters
    if isinstance(tiles, str):
        tiles = tiles.replace("'", "").split()

    # Get the burn data
    logger.info("Collecting MODIS burn data.")
    burn_data = BurnData(project_directory=project_directory, n_cores=n_cores)
    tiles = burn_data.get_burns(
        tiles=tiles,
        country=country,
        shape_file=shape_file,
        start_year=start_year,
        end_year=end_year
    )

    # Create a model builder object
    models = ModelBuilder(
        project_directory=project_directory,
        tiles=tiles,
        spatial_param=spatial_param,
        temporal_param=temporal_param,
        country=country,
        shape_file=shape_file,
        start_year=start_year,
        end_year=end_year,
        n_cores=n_cores
    )

    # Build event perimeters
    logger.info("Building event perimeter geometries.")
    gdf = models.build_events()

    # Sometimes no fire events are captured
    if gdf.shape[0] == 0:
        logger.info(
            f"No fire events found in tiles: {tiles} for years "
            f"{start_year} - {end_year}."
        ) # add new line to readme text file here: [Minor thought] Modis tile map: not all are covered by MCD64 (due to no fire) #128

    else:
        if land_cover_type==0: # add skip if the user didn't change the option
            msg = f"Skipping land classifciation: land_cover_type = {land_cover_type}"
            logger.info(msg)
        # Add land cover characteristic if requested
        else:
            lc_desc = ATTR_DESCS["landcover"][str(land_cover_type)]
            logger.info(f"Adding landcover attributes: {lc_desc}.")
            land_cover = LandCover(
                project_directory=project_directory,
                n_cores=n_cores
            )
            land_cover.get_land_cover(
                gdf=gdf,
                tiles=tiles,
                land_cover_type=land_cover_type
            )
            gdf = land_cover.add_land_cover_attributes(
                gdf=gdf,
                tiles=tiles,
                land_cover_type=land_cover_type
            )

        # Add ecoregion attributes
        ## Defaults are to not include eco_region
        if eco_region_type==None and eco_region_level==0:
            msg = f"Skipping ecoregion attributes (if this was unintended, check eco_region_type != {eco_region_type} and eco_region_level != {eco_region_level})"
            logger.info(msg)
        ## specifically account for requested types, if neither are used, throw a warning (or error?)
        ### I think add_eco_region_attributes can be consolidated (maybe removed) and we can directly call either the NA or WWF function - RE
        elif eco_region_type == "world":
            eco_desc = ATTR_DESCS["ecoregions"][eco_region_type]
            msg = f"Adding World Terrestrial Ecoregions from WWF."
            logger.info(msg)
            eco_region_data = EcoRegion(project_directory=project_directory)
            gdf = eco_region_data.add_eco_region_attributes(
                gdf=gdf,
                eco_region_type=eco_region_type
            )
        elif eco_region_type == "na" or eco_region_level>0:
            # add check to make sure eco region level is not 0. If it is, set it to 1, lowest level detail
            if eco_region_level==0:eco_region_level=1 
            eco_desc = ATTR_DESCS["ecoregions"][eco_region_type]
            msg = f"Adding NA ecoregions from North American ecoregions: {eco_desc}, Level {eco_region_level}."
            logger.info(msg)            
            eco_region_data = EcoRegion(project_directory=project_directory)
            eco_region_data.get_eco_region()
            gdf = eco_region_data.add_eco_region_attributes(
                gdf=gdf,
                eco_region_type=eco_region_type,
                eco_region_level=eco_region_level
            )
        else:
            logger.info("Eco region type was requested, but there seems to be a region problem") # this is logged if something other than na, world, or None is inputted.
            
        # Save the events
        models.save_data(
            gdf=gdf,
            project_name=project_name,
            project_directory=project_directory,
            start_year=start_year,
            end_year=end_year,
            daily=daily,
            shape_type=shape_type,
            csv_type=csv_type
        )

        # Done with processing, collect time and memory usage
        _, peak_memory = tracemalloc.get_traced_memory()
        peak_memory = round(peak_memory / 1024 ** 3, 2)
        tracemalloc.stop()
        end = time.perf_counter()
        seconds = end - start
        runtime = seconds / 60
        logger.info(f"Job completed in {runtime:.2f} minutes")
        logger.info(f"Peak memory usage: {peak_memory:.2f} GB")
        logger.info(f"Job completed in {runtime:.2f} minutes")
        logger.info(f"Peak memory usage: {peak_memory:.2f} GB")

        # Create a summary readme
        make_read_me(
            gdf=gdf,
            project_directory=project_directory,
            tiles=tiles,
            spatial_param=spatial_param,
            temporal_param=temporal_param,
            shapefile=shape_file,
            runtime=runtime,
            n_cores=n_cores,
            peak_memory=peak_memory,
            start_year=start_year,
            end_year=end_year,
            run_name=run_name,
            country=country,
        )

    # Remove intermediate files if requested
    
    if cleanup:
        logger.info(f"Cleaning up intermediate files in {project_directory}.")
        cleanup_intermediate_files(project_directory)

    return gdf
