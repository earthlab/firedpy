# -*- coding: utf-8 -*-
import argparse
import os
import shutil
import tempfile
import time
import warnings
import sys

from .functions import DataGetter, ModelBuilder
from .readme import makeReadMe

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
        event data frame. Spefify either "shp", "gpkg", or both. Shapefiles of both daily progression and overall
        event perimeters will be written to the "outputs/shapefiles" folder of
        the chosen project directory. These will be saved in the specified geopackage format
        (.gpkg), ERSI Shapefile format (.shp), or save them in both formats using the file basename of the fire event data frame (e.g.
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
    start_yr = ("""
        The first year of fired events.
        """)
    end_yr = ("""
        The last year of fired events.
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
    parser.add_argument("-shp_type", dest = "shp_type", help=shp_help, default = None)
    #parser.add_argument("--shpfile", action='store_true', help=shpf_help)
    parser.add_argument("-spatial", dest="spatial_param", default=5,
                        type=int, help=sp_help)
    parser.add_argument("-temporal", dest="temporal_param", default=11,
                        type=int, help=tmp_help)
    parser.add_argument("-aoi", "--names-list", nargs="+", dest="tiles",
                        default=["h08v04", "h09v04", "h10v04", "h11v04",
                                 "h12v04", "h13v04", "h08v05", "h09v05",
                                 "h10v05", "h11v05", "h13v04", "h08v05",
                                 "h09v05", "h10v05", "h11v05", "h12v05",
                                 "h08v06", "h09v06", "h10v06", "h11v06"],
                        help=tile_help)
    parser.add_argument("-daily", dest="daily", default="no", help=daily_help)
    parser.add_argument("-start_yr", dest="start_yr", type=int, default=None, help=start_yr)
    parser.add_argument("-end_yr", dest="end_yr", type=int, default=None, help=end_yr)


    if len(sys.argv) == 1:
        proj_dir = input("Enter project directory: ")
        if proj_dir == '':
            proj_dir=os.path.join(os.getcwd(), 'proj')
        if not os.path.exists(proj_dir):
            proj_dir=os.path.join(os.getcwd(), 'proj')
            os.makedirs(proj_dir)
        tilechoice = input("Would you like the fired product on a) continent, b) country, c) US state, or d)speific MODIS tiles?[a/b/c/d]: ")
        while True:
            if tilechoice not in ['a','b','c','d']:
                tilechoice = input("Enter a,b,c, or d: ")
            else:
                break
        if tilechoice == 'a':
            tilename = input("Please enter the continent name: /n")
            tiles = ["ref/continents/"+tilename +".gpkg", ]
            path = os.path.join("ref/continents/",tilename +".gpkg")
            while True:
                if os.path.isfile(path):
                    break
                else:
                    print("Not a valid choice. If you would like to see a list of avalible continents please visit: https://github.com/earthlab/firedpy/n")
                    tilename = input("Please enter the continent name: /n ")
                    tiles = ["ref/continents/"+tilename +".gpkg", ]
                    path = os.path.join("ref/continents/",tilename +".gpkg")

            if tilename ==  "north_america":
                ecoregion_type = 'na'
                ecoregion_level = 3
                landcover_type = 1
            else:
                ecoregion_type = 'world'
                landcover_type = None
                ecoregion_level = None

        if tilechoice == 'b':
            tilename = input("Please enter the country name:")
            tiles = ["ref/individual_countries/"+tilename +".gpkg", ]
            path = os.path.join("ref/individual_countries/",tilename +".gpkg")
            while True:
                if os.path.isfile(path):
                    break
                else:
                    print("Not a valid choice. If you would like to see a list of avalible countries please visit: https://github.com/earthlab/firedpy \n")
                    tilename = input("Please enter the country name: ")
                    tiles = ["ref/individual_countries/"+tilename +".gpkg", ]
                    path = os.path.join("ref/individual_countries/",tilename +".gpkg")


            na = ['united_States_of_america', 'canada', 'united_states_virgin_islands']
            if tilename in na:
                ecoregion_type = 'na'
                ecoregion_level = 3
                landcover_type = 1
            else:
                ecoregion_type = 'world'
                landcover_type = None
                ecoregion_level = None

        if tilechoice == 'c':
            tilename = input("Please enter the state name: ")
            tiles = ["ref/us_states/"+tilename +".gpkg", ]
            path = os.path.join("ref/us_states/",tilename +".gpkg")
            while True:
                if os.path.isfile(path):
                    break
                else:
                    print("Not a valid choice. If you would like to see a list of avalible states name please visit: https://github.com/earthlab/firedpy/n")
                    tilename = input("Please enter the state name: /n ")
                    tiles = ["ref/us_states/"+tilename +".gpkg", ]
                    path = os.path.join("ref/us_states/",tilename +".gpkg")

            ecoregion_type = 'na'
            ecoregion_level = 3
            landcover_type = 1

        if tilechoice == 'd':
            tiles = input("Please enter tiles as a list of characters (no quotes no spaces)(e.g., h08v04 h09v04 ...):")
            ecoregion_type = input("Please enter the ecoregion type [na or world]:")
            if ecoregion_type == 'na':
                ecoregion_level = 3
                landcover_type = 1
            else:
                ecoregion_level = None
                landcover_type =  None

        daily = input("Do you want to create the daily polygons or just the event-level perimeter for your analysis area (yes for daily/no for event-level): \n")
        if daily == '':
            daily = 'no'
        spatial_param = input("Please enter the number of cells (~463 m resolution) to search for neighboring burn detections. Defaults to 5 cells in all directions.\n")
        temporal_param = input("The number of days to search for neighboring burn detections. Defaults to 11 days between events.\n")
        if spatial_param == '':
            spatial_param = 5
        else:
            spatial_param = int(spatial_param)
        if temporal_param == '':
            temporal_param = 11
        else:
            temporal_param = int(temporal_param)

        shp_type = input("""Specify the format of the shapefile you want, [gpkg, shp, both, none] """)
        options = ["gpkg", "shp", "both", "none"]
        while True:
            if shp_type.lower() in options:
                break
            else:
                shp_type = input(""" Specified input was not valid \n Specify the format of the shapefile you want, "gpkg", "shp", "both", or "none" in all lower case """)
        if shp_type != 'none':
            shapefile = True
        elif shp_type == 'none':
            shapefile = False
            shp_type = None
        file_name = "fired_"+str(tilename)
        file_path = os.path.join(proj_dir,
                                     "outputs", "tables",
                                     file_name)
        start_yr = input("Enter the year you want to start or press enter for all dates: ")
        end_yr = input("Enter the year you want to end or press enter for all dates: ")
        if start_yr!='' and end_yr!='':
            start_yr = int(start_yr)
            end_yr = int(end_yr)
        else:
            start_yr = None
            end_yr = None
        input = 1
    else:

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
        shp_type = args.shp_type
        start_yr = args.start_yr
        end_yr = args.end_yr
        if shp_type != 'none':
            shapefile = True
        elif shp_type == 'none':
            shapefile = False
            shp_type = None

        input = 2
        name = str(tiles[0])
        nums=[str(0),str(1),str(4),str(3),str(4),str(5),str(6),str(7),str(8),str(9)]
        for i in nums:
            if i in name:
                input = 3
                tilename = tiles
        if input != 3:
            name = name.split("/")
            name = name[-1]
            name = name.split(".")
            tilename = name[0]

            # Assign the temporary file name including the spatial and temporal parameters
            file_path = os.path.join(args.proj_dir,
                                         "outputs", "tables",
                                         args.file_name+"_"+str(tilename))
        else:
            # Assign the temporary file name including the spatial and temporal parameters
            file_path = os.path.join(args.proj_dir,
                                         "outputs", "tables",
                                         args.file_name)



    # Transfer the lookup tables
    if landcover_type:
        lookup = os.path.join(os.getcwd(), 'ref', 'landcover',
                              'MCD12Q1_LegendDesc_Type{}.csv'.format(str(landcover_type)))
                              
        new_path = os.path.join(proj_dir, 'tables', 'landcover')
        if not os.path.exists(new_path):
            os.makedirs(new_path)
        new_file = os.path.join(new_path, 'MCD12Q1_LegendDesc_Type{}.csv'.format(str(landcover_type)))
        shutil.copy(lookup, new_file)
                              
                              
                              
     # Get ecoregions if requested, use local file first
    if ecoregion_type or ecoregion_level:

        new_path = os.path.join(proj_dir, 'shapefiles', 'ecoregion')

        if not os.path.exists(new_path):
            os.makedirs(new_path)

        if ecoregion_type == 'world':
            fname = 'wwf_terr_ecos.gpkg'
            lookup = os.path.join(os.getcwd(), 'ref', 'world_ecoregions', fname)
        elif ecoregion_type == 'na' or ecoregion_level:
            fname = 'NA_CEC_Eco_Level3.gpkg'
            lookup = os.path.join(os.getcwd(), 'ref', 'us_eco', fname)
            print(lookup)
        try:
            new_file = os.path.join(new_path, fname)
            shutil.copy(lookup, new_file)
        except Exception:
            data.getEcoregion(ecoregion_level)
                              

    # Make sure the project directory exists
    if not os.path.exists(proj_dir):
        os.makedirs(proj_dir)

    # Create data object
    data = DataGetter(proj_dir, start_yr, end_yr)

    # Assign target MODIS tiles to the data object
    if os.path.splitext(tiles[0])[1] in [".shp", ".gpkg"]:
        shp = tiles[0]
        print("Filtering for MODIS tiles that intersect \n    " + shp)
        data.shapeToTiles(shp)
        tiles = data.tiles
    else:
        data.tiles = tiles
        shp = tiles[0]

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


    # Grab the basename for the output file name
    file_base = os.path.basename(file_path)

    # Add date range to the file names before exporting final data frame
    date_range = []
    for root, dirs, files in os.walk(os.path.join(proj_dir, 'rasters', 'burn_area', 'hdfs')):
        for f in files:
            if "MCD64A1" in f:
                dr = int(f.split('.')[1][1:])
                date_range.append(dr)
    last_date = sorted(date_range)[-1]
    first_date = sorted(date_range)[1]
    file_name = os.path.join(os.path.dirname(file_path), file_base+"_to"+str(last_date))



    # Create Model Builder object
    models = ModelBuilder(file_name=file_name,
                          proj_dir=proj_dir,
                          tiles=tiles,
                          shp=shp,
                          spatial_param=spatial_param,
                          temporal_param=temporal_param,
                          landcover_type=landcover_type,
                          ecoregion_type=ecoregion_type,
                          ecoregion_level=ecoregion_level,
                          daily=daily,
                          shapefile=shapefile,
                          shp_type=shp_type)

    # Now go ahead and create the events (Memory's a bit tight for parallel)
    models.buildEvents()

    # Now add fire attributes to this table
    models.buildFireAttributes()

    # And now build the polygons
    pref = file_base+"_to"+str(last_date)
    daily_shp_file = "_".join([pref, "daily"])
    event_shp_file = "_".join([pref, "events"])

    if shp_type == "shp":
        daily_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
                                      daily_shp_file + ".shp")
        event_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
                                      event_shp_file + ".shp")
    if shp_type == "both":
        daily_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
                                      daily_shp_file + ".shp")
        daily_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      daily_shp_file + ".gpkg")
        event_shp_path_shp = os.path.join(proj_dir, "outputs", "shapefiles",
                                      event_shp_file + ".shp")
        event_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      event_shp_file + ".gpkg")
    else:
        daily_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      daily_shp_file + ".gpkg")
        event_shp_path = os.path.join(proj_dir, "outputs", "shapefiles",
                                      event_shp_file + ".gpkg")
        daily_shp_path_shp = ''
        event_shp_path_shp = ''


    models.buildPolygons(daily_shp_path=daily_shp_path,
                         event_shp_path=event_shp_path,
                         daily_shp_path_shp = daily_shp_path_shp ,
                         event_shp_path_shp =daily_shp_path_shp )
    makeReadMe(proj_dir, tilename, file_base, input, first_date, last_date, ecoregion_type, ecoregion_level, landcover_type, daily, spatial_param, temporal_param, shapefile, shp_type)
    # Print the time it took
    end = time.perf_counter()
    seconds = end - start
    minutes = seconds/60
    print("Job completed in {} minutes".format(round(minutes, 2)))



if __name__ == "__main__":
    main()
