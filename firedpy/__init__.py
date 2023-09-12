name = "firedpy"
DATA_HELP = ("""
    The project directory you would like to use for  and output
    data files. Defaults to a temporary directory 'firedpy/proj'.
    """)
FILE_HELP = ("""
    The file name of the resulting dataframe. This will be saved in
    the "outputs/tables" folder of the chosen project directory. Defaults
    to "fired_events.csv" and "fired_daily.csv" if daily data is requested.
    """)
DAILY_HELP = ("""
    You may specify whether to create the daily polygons or just the event-level perimeter
    for your analysis area. Options are "yes" (to create the daily polygons and the event polygons),
    "no" (create the event level only).
    """)
ECO_HELP = ("""
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
LC_HELP = ("""
    To include land cover as an attribute, provide a number corresponding
    with a MODIS/Terra+Aqua Land Cover (MCD12Q1) category followed with :username:password of your NASA's Earthdata service
    account. Available land cover categories:
        1: IGBP global vegetation classification scheme,
        2: University of Maryland (UMD) scheme,
        3: MODIS-derived LAI/fPAR scheme.

    If you do not have an account register at https://urs.earthdata.nasa.gov/home. Defaults to none.
    """)
SHP_HELP = ("""
    Provide this option if you would like to build shapefiles from the
    event data frame. Spefify either "shp", "gpkg", or both. Shapefiles of both daily progression and overall
    event perimeters will be written to the "outputs/shapefiles" folder of
    the chosen project directory. These will be saved in the specified geopackage format
    (.gpkg), ERSI Shapefile format (.shp), or save them in both formats using the file basename of the fire event data frame (e.g.
    'modis_events_daily.gpkg' and 'modis_events.gpkg')
    """)
SP_HELP = ("""
    The number of cells (~463 m resolution) to search for neighboring burn
    detections. Defaults to 5 cells in all directions.
    """)
TILE_HELP = ("""
    You may specify the tiles as a list of characters (no quotes no spaces)
    (e.g., h08v04 h09v04 ...) or leave this blank to default to tiles
    covering the Contiguous United States. Specify "all" to use all
    available MODIS tiles. Alternatively, provide a path to a shapefile
    with either a ".shp" or ".gpkg" extension to use intersecting MODIS
    tiles. In the firedpy directory, you can access any of the 50 states by specifying ref/us_states/state_name.gpkg,
    all 7 continents by specifying ref/continents/continent_name.gpkg, and a country by ref/individual_countries/country_name.gpkg.
    A list of all avalible countries can be found in the ReadMe.
    """)
TMP_HELP = ("""
    The number of days to search for neighboring burn detections. Defaults
    to 11 days between events.
    """)
START_YR = ("""
    The first year of fired events.
    """)
END_YR = ("""
    The last year of fired events.
    """)
FULL_CSV = ("""
    If included full attribute table will exported to csv. If not included only x and y coordinates, event date, and
    event id will be exported to a csv.""")