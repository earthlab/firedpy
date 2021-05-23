# -*- coding: utf-8 -*-
import argparse
import os
import shutil
import tempfile
import time
import warnings

from .functions import DataGetter, ModelBuilder

warnings.filterwarnings("ignore", category=FutureWarning)


def main():
    # Start the timer (seconds)
    start = time.perf_counter()

    # Call help statements
    data_help = ("""
        The project directory you would like to use for input and output
        data files. Defaults to a temporary directory 'firedpy/proj'.
        """)
    file_help = ("""
        The file name of the resulting dataframe. This will be saved in
        the "outputs/tables" folder of the chosen project directory. Defaults
        to "fired_events.csv" and "fired_daily.csv" if daily data is requested.
        """)
    daily_help = ("""
        You may specify whether to create the daily polygons or just the event-level perimeter
        for your analysis area. Options are "yes" (to create the daily polygons and the event polygons),
        "no" (create the event level only).
        """)
    eco_help = ("""
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
    lc_help = ("""
        To include land cover as an attribute, provide a number corresponding
        with a MODIS/Terra+Aqua Land Cover (MCD12Q1) category. To do so you
        will have to register at NASA's Earthdata service
        (https://urs.earthdata.nasa.gov/home) and enter your user name and
        password when prompted. Available land cover categories:
            1: IGBP global vegetation classification scheme,
            2: University of Maryland (UMD) scheme,
            3: MODIS-derived LAI/fPAR scheme,
            4: MODIS-derived Net Primary Production (NPP) scheme,
            5: Plant Functional Type (PFT) scheme.
        Defaults to none.
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
        The number of cells (~463 m resolution) to search for neighboring burn
        detections. Defaults to 5 cells in all directions.
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
                        default=os.path.join(os.getcwd(), 'proj'), help=data_help)
    parser.add_argument("-file_name", dest="file_name",
                        default="fired", help=file_help)
    parser.add_argument("-ecoregion_type", dest="ecoregion_type", default=None,
                        help=eco_help)
    parser.add_argument("-ecoregion_level", dest="ecoregion_level", type=int,
                        default=None, help=eco_help)
    parser.add_argument("-landcover_type", dest="landcover_type", type=int,
                        default=None, help=lc_help)
    parser.add_argument("--shapefile", action='store_true', help=shp_help)
    parser.add_argument("-spatial", dest="spatial_param", default=5,
                        type=int, help=sp_help)
    parser.add_argument("-temporal", dest="temporal_param", default=11,
                        type=int, help=tmp_help)
    parser.add_argument("-tiles", "--names-list", nargs="+", dest="tiles",
                        default=["h08v04", "h09v04", "h10v04", "h11v04",
                                 "h12v04", "h13v04", "h08v05", "h09v05",
                                 "h10v05", "h11v05", "h13v04", "h08v05",
                                 "h09v05", "h10v05", "h11v05", "h12v05",
                                 "h08v06", "h09v06", "h10v06", "h11v06"],
                        help=tile_help)
    parser.add_argument("-daily", dest="daily", default="no", help=daily_help)

    # Parse argument responses
    args = parser.parse_args()
    proj_dir = args.proj_dir
    ecoregion_type = args.ecoregion_type
    ecoregion_level = args.ecoregion_level
    landcover_type = args.landcover_type
    daily = args.daily
    spatial_param = args.spatial_param
    temporal_param = args.temporal_param
    tiles = args.tiles
    shapefile = args.shapefile

    # Assign the temporary file name including the spatial and temporal parameters
    file_name = os.path.join(args.proj_dir,
                             "outputs", "tables",
                             args.file_name + "_events.csv")

    # Transfer the lookup tables
    if landcover_type:
        lookup = os.path.join(os.getcwd(), 'ref', 'landcover',
                              'MCD12Q1_LegendDesc_Type{}.csv'.format(str(landcover_type)))
        new_path = os.path.join(proj_dir, 'tables', 'landcover')
        if not os.path.exists(new_path):
            os.makedirs(new_path)
        new_file = os.path.join(new_path, 'MCD12Q1_LegendDesc_Type{}.csv'.format(str(landcover_type)))
        shutil.copy(lookup, new_file)

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

    # Get ecoregions if requested
    if ecoregion_type or ecoregion_level:
        data.getEcoregion(ecoregion_level)

    # Add date range to the file names before exporting final data frame
    date_range = []
    for root, dirs, files in os.walk(os.path.join(proj_dir, "rasters", "burn_area", "hdfs")):
        for f in files:
            dr = int(f.split('.')[1][1:])
            date_range.append(dr)
    last_date = sorted(date_range)[-1]
    file_name = file_name[:-4]+"_"+str(last_date)+".csv"

    # Create Model Builder object
    models = ModelBuilder(file_name=file_name,
                          proj_dir=proj_dir,
                          tiles=tiles,
                          spatial_param=spatial_param,
                          temporal_param=temporal_param,
                          landcover_type=landcover_type,
                          ecoregion_type=ecoregion_type,
                          ecoregion_level=ecoregion_level,
                          daily=daily,
                          shapefile=shapefile)

    # Now go ahead and create the events (Memory's a bit tight for parallel)
    models.buildEvents()

    # Now add fire attributes to this table
    models.buildFireAttributes()

    # And now build the polygons
    file_base = os.path.basename(args.file_name)

    daily_shp_file = "_".join([file_base, "daily"])
    event_shp_file = "_".join([file_base, "events"])

    daily_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                  daily_shp_file + ".gpkg")
    event_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                  event_shp_file + ".gpkg")

    models.buildPolygons(daily_shp_path=daily_shp_path,
                         event_shp_path=event_shp_path)

    # Print the time it took
    end = time.perf_counter()
    seconds = end - start
    minutes = seconds/60
    print("Job completed in {} minutes".format(round(minutes, 2)))


if __name__ == "__main__":
    main()
