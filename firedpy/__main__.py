# -*- coding: utf-8 -*-
import argparse
import os
import tempfile
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
from functions import buildEvents, Data_Getter

def main():
    # Call help statements
    tile_help = ("""
        You may specify the tiles as a list of characters (no quotes no spaces)
        (e.g., h08v04 h09v04 ...) or leave this blank to default to tiles
        covering the Contiguous United States. Specify 'all' to use all
        available MODIS tiles.
        """)
    shp_help = ("""
        You may also specify a path to a shapefile to set the MODIS grid tiles
        to use. Use local shapefiles with the .shp extension or a url to a zip
        file of shapefile elements. Defaults to None.
        """)
    sp_help = ("""
        The number of cells (463 m2 each) to search for neighboring burn
        detections. Defaults to 5 cells in all directions.
        """)
    tmp_help = ("""
        The number of days to search for neighboring burn detections. Defaults
        to 11 days between events.
        """)
    dest_help = ("""
        Specify the location and filename of the resulting dataframe. Defaults
        to a file called modis_events.csv in a temporary directory.
        """)
    data_help = ("""
        specify the directory you'd like to use for the required shape and
        and raster data files. Defaults to a temporary directory.
        """)

    # Provide arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-dest', dest='dest',
                        default=os.path.join(tempfile.mkdtemp(),
                                             'modis_events.csv'),
                        help=dest_help)
    parser.add_argument('-tiles', '--names-list', nargs="+", dest='tiles',
                        default=["h08v04", "h09v04", "h10v04", "h11v04",
                                 "h12v04", "h13v04", "h08v05", "h09v05",
                                 "h10v05", "h11v05", "h13v04", "h08v05",
                                 "h09v05", "h10v05", "h11v05", "h12v05",
                                 "h08v06", "h09v06", "h10v06", "h11v06"],
                        help=tile_help)
    parser.add_argument('-shapefile', dest='shp', default=None, help=shp_help)
    parser.add_argument('-spatial_param', dest='spatial_param', default=5,
                        type=int, help=sp_help)
    parser.add_argument('-temporal_param', dest='temporal_param', default=11,
                        type=int, help=tmp_help)
    parser.add_argument('-data_dir', dest='data_dir',
                        default=tempfile.mkdtemp(), help=data_help)

    # Parse argument responses
    args = parser.parse_args()
    tiles = args.tiles
    spatial_param = args.spatial_param
    temporal_param = args.temporal_param
    args = parser.parse_args()
    tiles = args.tiles
    shp = args.shp
    dest = args.dest
    data_dir = args.data_dir

    # Create data object (you can specify which tiles as an attribute)
    data = Data_Getter(data_dir)
    if shp:
        print("Filtering for MODIS tiles that intersect \n    " + shp)
        data.shapeToTiles(shp)
        tiles = data.tiles
    else:
        data.tiles = tiles

    # Get all of the MODIS burn area hdfs (this needs to be parallelized)
    data.getBurns()

    # Now go ahead and create the events (Memory's a bit tight for parallel)
    buildEvents(dest=dest, data_dir=data_dir, tiles=tiles,
                spatial_param=spatial_param, temporal_param=temporal_param)

    # And build the polygons
    # ...

if __name__ == '__main__':
    main()