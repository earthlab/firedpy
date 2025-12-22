HELP_TEXT = {
    "name": "src",
    "cleanup": (
        """
        If set then the burn area and landcover files will be removed after
        each run to save disk space in between multiple runs.
        """
    ),
    "daily": (
        """
        You may specify whether to create the daily polygons or just the
        event-level perimeter for your analysis area. If this flag is set, the
        daily and event polygons will be created, otherwise only the event
        level.
        """
    ),
    "data": (
        """
        The project directory you would like to use for  and output data files.
        Defaults to a temporary directory 'firedpy/data/outputs'.
        """
    ),
    "file": (
        """
        The file name of the resulting dataframe. This will be saved in the
        'outputs/tables' folder of the chosen project directory. Defaults to
        'fired_events.csv' and 'fired_daily.csv' if daily data is requested.
        """
    ),
    "earthdata": (
        """
        Please input your NASA Earthdata username and password in order to
        download the land cover data. If you do not have an Earthdata account,
        you can register at https://urs.earthdata.nasa.gov/. To avoid seeing
        this prompt again, you can set the FIREDPY_ED_USER and FIREDPY_ED_PWD
        environment variables.
        """
    ),
    "eco": (
        """
        You can specify the ecoregion type as either 'world' or 'na':\n
        'world' = World Terrestrial Ecoregions (World Wildlife Fund (WWF))
        'na' = North American ecoregions (Omernick, 1987)

        Most common (modal) ecoregion across the event is used.

        Further, to associate each event with North American ecoregions
        (Omernick, 1987) you may provide a number corresponding to an ecoregion
        level. Ecoregions are retrieved from www.epa.gov and levels I through
        IV are available. Levels I and II were developed by the North American
        Commission for Environmental Cooperation. Levels III and IV were
        developed by the United States Environmental Protection Agency. For
        events with more than one ecoregion, the most common value will be
        used. Defaults to None.
        """
    ),
    "full_csv": (
        """
        If set full attribute table will exported to csv. If not included only
        x and y coordinates, event date, and event id will be exported to a
        csv.
        """
    ),
    "lc": (
        """
        To include land cover as an attribute, provide a number corresponding
        with a MODIS/Terra+Aqua Land Cover (MCD12Q1) category followed with
        :username:password of your NASA's Earthdata service account. Available
        land cover categories:\n

            1: IGBP global vegetation classification scheme,
            2: University of Maryland (UMD) scheme,
            3: MODIS-derived LAI/fPAR scheme.

        If you do not have an account register at
        https://urs.earthdata.nasa.gov/home. Defaults to None.
        """
    ),
    "n_cores": "Number of cores to use for parallel processing.",
    "shp": (
        """
        Provide this option if you would like to build shapefiles from the
        event data frame. Specify either "shp", "gpkg", or both. Shapefiles of
        both daily progression and overall event perimeters will be written to
        the 'outputs/shapefiles' folder of the chosen project directory. These
        will be saved in the specified geopackage format (.gpkg), ERSI
        Shapefile format (.shp), or save them in both formats using the file
        basename of the fire event data frame (e.g. 'modis_events_daily.gpkg'
        and 'modis_events.gpkg').
        """
    ),
    "tile": (
        """
        You may specify the tiles as a list of characters (no quotes no
        spaces) (e.g., 'h08v04', 'h09v04', etc.). Alternatively, provide a path
        to a shapefile with either a '.shp' or '.gpkg' extension to use its
        intersecting MODIS tiles. In the `src` directory, you can access any of
        the 50 states by specifying 'ref/us_states/state_name.gpkg', all 7
        continents by specifying 'ref/continents/continent_name.gpkg', and a
        country by 'ref/individual_countries/country_name.gpkg'. A list of all
        available countries can be found in the ReadMe.
        """
    ),
    "sp": (
        """
        The number of cells (~463 m resolution) to search for neighboring burn
        detections. Defaults to 5 cells in all directions.
        """
    ),
    "tile_name": (
        """
        The name of the tile you would like to choose based on your tile
        choice.
        """
    ),
    "tmp": (
        """
        The number of days to search for neighboring burn detections. Defaults
        to 11 days between events.
        """
    ),
    "end_year": "The last year of fired events.",
    "start_year": "The first year of fired events.",
}
