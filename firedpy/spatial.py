import os

from typing import List

import geopandas as gpd
import pandas as pd

from firedpy import DATA_DIR


def shape_to_tiles(shape_path: str) -> List[str]:
    """
    Set or reset the tile list using a shapefile. Where shapes intersect
    with the modis sinusoidal grid determines which tiles to use.
    """

    # MODIS CRS retrieved from a single HDF file
    out_crs = '''PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",
    DATUM["Not specified (based on custom spheroid)",
    SPHEROID["Custom spheroid",6371007.181,0]],
    PRIMEM["Greenwich",0],
    UNIT["degree",0.0174532925199433,
    AUTHORITY["EPSG","9122"]]],
    PROJECTION["Sinusoidal"],
    PARAMETER["longitude_of_center",0],
    PARAMETER["false_easting",0],
    PARAMETER["false_northing",0],
    UNIT["Meter",1],
    AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''

    grid_path = DATA_DIR.joinpath("modis_grid.gpkg")
    modis_grid = gpd.read_file(grid_path)
    modis_grid.set_crs(out_crs, inplace=True, allow_override=True)

    # Read in the input shapefile and reproject to MODIS sinusoidal
    source = gpd.read_file(shape_path)
    source.to_crs(out_crs, inplace=True)

    # Left join shapefiles with source shape as the left
    shared = gpd.sjoin(source, modis_grid, how="left").dropna()
    shared["h"] = shared["h"].apply(lambda x: "h{:02d}".format(int(x)))
    shared["v"] = shared["v"].apply(lambda x: "v{:02d}".format(int(x)))
    shared["tile"] = shared["h"] + shared["v"]
    tiles = pd.unique(shared["tile"].values)

    return tiles
