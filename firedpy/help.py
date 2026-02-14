ATTR_HELP = {
    "eco_region_level": {
        1: "Level I - Least Detailed",
        2: "Level III - Mid-level Detail",
        3: "Level III - Most Detailed"
    },
    "eco_region_type": {
        "world": "World Terrestrial Ecoregions (World Wildlife Fund)",
        "na": "North American Ecoregions (Omernick, 1987)"
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
        "gpkg": "GeoPackage",
        "both": "Both ESRI Shapefile and GeoPackage"
    }
}


CLI_HELP = {
    "cleanup": (
        """
        Cleanup. If set then the burn area and landcover files will be removed
        after each run to save disk space in between multiple runs.
        """
    ),
    "country": (
        """
        Country. The name of a country to use as a study area. If not provided,
        either a `tiles` or `shape_file` parameter is required. Country
        overrides `tiles` parameter.
        """
    ),
    "daily": (
        """
        Daily. Create daily polygons in addition to event-level perimeters for
        your analysis area.
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
        values and confirm selections before proceeding. Defaults to False.
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
        """
    ),
    "n_cores": (
        """
        N Cores. Number of cores to use for parallel processing. A value of 0
        will use all available cores.
        """
    ),
    "project_directory": (
        """
        Project directory path. Inputs and outputs will all be written here.
        Required (Use "." for the present directory).
        """
    ),
    "project_name": (
        """
        A name used to identify the output files of this project. Defaults to
        None, which will use the name of the parent run directory.
        """
    ),
    "shape_file": (
        """
        Shapefile. A shapefile representing the target study area. If not
        provided, either a `tiles` or `country` parameter is required. This
        overrides the `country` and `tiles` parameters. Defaults to None.
        """
    ),
    "shape_type": (
        """
        Shapefile Type. The file format for the output event GeodataDataFrame.
        Specify "shp" for an ESRI Shapefile, "gpkg" for a Geopackage, or "both"
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
        Tiles. A string representing a single MODIS tile (e.g., 'h08v04') or a
        string representing multiple tiles separated by spaces (e.g.,
        'h08v04 h09v04'). If not provided, a `country` or `shape_file`
        parameter is required. Defaults to None
        """
    )
}
