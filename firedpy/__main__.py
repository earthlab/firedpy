# -*- coding: utf-8 -*-
import argparse
from .functions import DataGetter, ModelBuilder
import os
import tempfile
import time
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# Checkout http://effbot.org/zone/stupid-exceptions-keyboardinterrupt.htm

def main():
    # Start the timer (seconds)
    start = time.perf_counter()

    # Call help statements
    data_help = ("""
        The project directory you would like to use for input and output
        data files. Defaults to a temporary directory.
        """)
    dest_help = ("""
        The filename of the resulting dataframe. This will be saved in
        the "outputs/tables" folder of the chosen project directory. Defaults
        to "modis_events.csv".
        """)
    eco_help = ("""
        Provide this option to associate each event with an ecoregion. 
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
        The number of cells (463 m2 each) to search for neighboring burn
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
                        default=tempfile.mkdtemp(), help=data_help)
    parser.add_argument("-dest", dest="dest",
                        default="modis_events.csv",
                        help=dest_help)
    parser.add_argument("-ecoregion", dest="ecoregion", default=None,
                        help=lc_help)
    parser.add_argument("-landcover_type", dest="landcover_type", default=None,
                        help=lc_help)
    parser.add_argument("--shapefile", action='store_true', help=shp_help)
    parser.add_argument("-spatial_param", dest="spatial_param", default=5,
                        type=int, help=sp_help)
    parser.add_argument("-temporal_param", dest="temporal_param", default=11,
                        type=int, help=tmp_help)
    parser.add_argument("-tiles", "--names-list", nargs="+", dest="tiles",
                        default=["h08v04", "h09v04", "h10v04", "h11v04",
                                 "h12v04", "h13v04", "h08v05", "h09v05",
                                 "h10v05", "h11v05", "h13v04", "h08v05",
                                 "h09v05", "h10v05", "h11v05", "h12v05",
                                 "h08v06", "h09v06", "h10v06", "h11v06"],
                        help=tile_help)

    # Parse argument responses
    args = parser.parse_args()
    proj_dir = args.proj_dir
    ecoregion = args.ecoregion
    landcover_type = args.landcover_type
    dest = os.path.join(proj_dir, "outputs", "tables", args.dest)
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

    # Get all of the MODIS burn area hdfs
    try:
        data.getBurns()
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        template = "\nDownload failed: error type {0}:\n{1!r}"
        message = template.format(type(e).__name__, e.args)
        print(message)

    # Get land cover if requested
    if landcover_type:
        data.getLandcover(landcover_type)

    # Get ecoregions if requested
    if ecoregion:
        data.getEcoregion(ecoregion)

    # Create Model Builder object
    models = ModelBuilder(dest=dest,
                          proj_dir=proj_dir,
                          tiles=tiles,
                          spatial_param=spatial_param,
                          temporal_param=temporal_param,
                          landcover_type=landcover_type,
                          ecoregion=ecoregion)

    # Now go ahead and create the events (Memory's a bit tight for parallel)
    models.buildEvents()

    # Now add attributes to this table
    lc_dir = data.landcover_mosaic_path
    eco_dir = None
    models.buildAttributes(lc_dir=lc_dir, eco_dir=eco_dir)

    # And build the polygons
    if shapefile:
        # Use arguments for shapefile source and destination file paths
        file_base = os.path.splitext(os.path.basename(dest))[0]
        daily_shp_file = "_".join([file_base, "daily"])
        daily_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      daily_shp_file + ".gpkg")
        event_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      file_base + ".gpkg")
        models.buildPolygons(daily_shp_path=daily_shp_path,
                             event_shp_path=event_shp_path)

    # Print the time it took
    end = time.perf_counter()
    seconds = end - start
    minutes = seconds/60
    print("Job completed in {} minutes".format(round(minutes, 2)))

if __name__ == "__main__":
    main()
