# -*- coding: utf-8 -*-
import argparse
from .functions import buildEvents, buildPolygons, DataGetter
import os
import tempfile
import time
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

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
    dest = os.path.join(proj_dir, "output", "tables", args.dest)
    spatial_param = args.spatial_param
    temporal_param = args.temporal_param
    tiles = args.tiles
    shapefile = args.shapefile

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

    # Get all of the MODIS burn area hdfs (this needs to be parallelized)
    data.getBurns()

    # Now go ahead and create the events (Memory"s a bit tight for parallel)
    buildEvents(dest=dest, data_dir=proj_dir, tiles=tiles,
                spatial_param=spatial_param, temporal_param=temporal_param)

    # And build the polygons
    if shapefile:
        # Use arguments for shapefile source and destination file paths
        file_base = os.path.splitext(os.path.basename(dest))[0]
        daily_shp_file = "_".join([file_base, "daily"])
        daily_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      daily_shp_file + ".gpkg")
        event_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      file_base + ".gpkg")
        shp_src = dest
        buildPolygons(src=shp_src, daily_shp_path=daily_shp_path,
                      event_shp_path=event_shp_path, data_dir=proj_dir)

    # Print the time it took
    end = time.perf_counter()
    seconds = end - start
    minutes = seconds/60
    print("Job completed in {} minutes".format(round(minutes, 2)))

if __name__ == "__main__":
    main()
