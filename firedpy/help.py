ATTR_HELP = {
    "eco_region_level": {
        1: "Level I - Least Detailed",
        2: "Level III - Mid-level Detail",
        3: "Level III - Most Detailed"
    },
    "eco_region_type": {
        "world": "World Terrestrial Ecoregions (World Wildlife Fund)",
        "na": "North American ecoregions (Omernick, 1987)"
    },
    "land_cover_type": {
        1: "International Geosphere-Biosphere Programme (IGBP) scheme",
        2: "University of Maryland (UMD) scheme",
        3: "MODIS-derived Leaf Area Index (LAI/fPAR) scheme",
        4: "MODIS-derived Net Primary Production (NPP) scheme",
        5: "Plant Functional Type (PFT) scheme"
    },
    "shape_type": {
        "shp": "ESRI Shapefile",
        "gpkg": "Geopackage",
        "both": "Both ESRI Shapefile and Geopackage"
    }
}


CLI_HELP = {
    "cleanup": (
        """
        Cleanup. If set then the burn area and landcover files will be removed
        after each run to save disk space in between multiple runs.
        """
    ),
    "daily": (
        """
        Daily. Create daily polygons in addition to event-level perimeters for
        your analysis area.
        """
    ),
    "data": (
        """
        Data. The project directory you would like to use for  and output data
        files. Defaults to a temporary directory 'firedpy/data/outputs'.
        """
    ),
    "country": (
        """
        The name of a country to use as a study area. Defaults to None.
        """
    ),
    "eco_region_type": (
        """
        Ecoregion Type. Specify the ecoregion type as either 'world' or 'na':\n
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
    "eco_region_level": (
        """
        Ecoregion Level. The desired Ecoregions level from the North American
        Commission for Environmental Cooperation (CEC). Levels 1 to 3 are
        available, with level 1 representing the broadest scale and level III
        representing the most detailed. Defaults to 1.
        """
    ),
    "end_year": "Last Year. The last year of fired events.",
    "file": (
        """
        File Name. The file name of the resulting dataframe. This will be saved
        in the 'outputs/tables' folder of the chosen project directory.
        Defaults to 'fired_events.csv' and 'fired_daily.csv' if daily data is
        requested.
        """
    ),
    "full_csv": (
        """
        Full CSV. Include a full set of attributes in the output CSV. If not
        included only x and y coordinates, event date, and event id will be
        exported to a csv.
        """
    ),
    "interactive": (
        """
        Interactive Mode. Firedpy will prompt the user for all parameter
        values. Defaults to False.
        """
    ),
    "land_cover_type": (
        """
        Land Cover Type. The target land cover source used to characterize each
        fire event. Provide a number corresponding with a MODIS/Terra+Aqua
        Land Cover (MCD12Q1). Available land cover categories:\n

            1: IGBP global vegetation classification scheme\n
            2: University of Maryland (UMD) scheme\n
            3: MODIS-derived LAI/fPAR scheme\n
            4: Annual BIOME-Biogeochemical Cycles (BGC)\n
            5: Annual Plant Functional Types (PFT)\n

        An Earthdata account is required. If you don't have one, register at
        https://urs.earthdata.nasa.gov/home. If you have an account, firedpy
        will use any of the authentication options outlined here:

        https://earthaccess.readthedocs.io/en/latest/user_guide/authenticate/

        Defaults to None.
        """
    ),
    "n_cores": (
        """
        N Cores. Number of cores to use for parallel processing. Defaults to
        0 or all cores.
        """
    ),
    "project_directory": (
        """
        Project directory path. Inputs and outputs will all be written here.
        Required (Use "." for the present directory).
        """
    ),
    "shape_file": (
        """
        Shapefiles. Provide this option if you would like to build shapefiles
        from the event data frame. Specify either "shp", "gpkg", or both.
        Shapefiles of both daily progression and overall event perimeters will
        be written to the 'outputs/shapefiles' folder of the chosen project
        directory. These will be saved in the specified geopackage format
        (.gpkg), ERSI Shapefile format (.shp), or save them in both formats
        using the file basename of the fire event data frame (e.g.
        'modis_events_daily.gpkg' and 'modis_events.gpkg').
        """
    ),
    "shape_type": (
        """
        Shapefile Type. The file format for the output event GeodataDataFrame.
        Specify "shp" for an ESRI Shapefil, "gpkg" for a Geopackage, or "both"
        to write both formats to file. Shapefiles of both daily
        progression and overall event perimeters will be written to the
        'outputs/shapefiles' folder of the chosen project directory.
        """
    ),
    "start_year": "First Year. The first year of fired events.",
    "spatial_param": (
        """
        Spatial Parameter. The number of cells (~463 m resolution) to search
        for neighboring burn detections. Defaults to 8 cells in all directions.
        """
    ),
    "temporal_param": (
        """
        Temporal Parameter. The number of days to search for neighboring burn
        detections. Defaults to 3 days between events.
        """
    ),
    "tile_name": (  #<--------------------------------------------------------- This doesn't make sense
        """
        Tile Name. The name of the tile you would like to choose based on your
        tile choice.
        """
    ),
    "tiles": (
        """
        Tiles. A list of MODIS tile IDs as a list of characters (no quotes no
        spaces) (e.g., ['h08v04', 'h09v04']). Alternatively, provide a path
        to a shapefile with either a '.shp' or '.gpkg' extension to use its
        intersecting MODIS tiles. In the `src` directory, you can access any of
        the 50 states by specifying 'ref/us_states/state_name.gpkg', all 7
        continents by specifying 'ref/continents/continent_name.gpkg', and a
        country by 'ref/individual_countries/country_name.gpkg'. A list of all
        available countries can be found in the ReadMe.
        """
    )
}
