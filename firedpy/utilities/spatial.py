"""Spatial data manipulation utilities."""
import os

import geopandas as gpd
import pandas as pd

from osgeo import gdal
from firedpy import DATA_DIR


# MODIS CRS retrieved from a single HDF file
MODIS_CRS = (
    """
    PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",
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
    AXIS["Easting",EAST],AXIS["Northing",NORTH]]
    """
)


def shape_to_tiles(shape_path):
    """Set or reset the tile list using a shapefile.

    NOTE: Where shapes intersect with the modis sinusoidal grid determines
    which tiles to use.

    Paramters
    ---------
    shape_path : str
        File path to the target shapefile.

    Returns
    -------
    list[str] : A list of strings representing MODIS tiles that intersect with
        with the the target shapefile.
    """
    # Read in the MODIS grid
    grid_path = DATA_DIR.joinpath("modis_grid.gpkg")
    modis_grid = gpd.read_file(grid_path)
    modis_grid.set_crs(MODIS_CRS, inplace=True, allow_override=True)

    # Read in the input shapefile and reproject to MODIS sinusoidal
    source = gpd.read_file(shape_path)
    source.to_crs(MODIS_CRS, inplace=True)

    # Left join shapefiles with source shape as the left
    shared = gpd.sjoin(source, modis_grid, how="left").dropna()
    shared["h"] = shared["h"].apply(lambda x: "h{:02d}".format(int(x)))
    shared["v"] = shared["v"].apply(lambda x: "v{:02d}".format(int(x)))
    shared["tile"] = shared["h"] + shared["v"]
    tiles = pd.unique(shared["tile"].values)

    return tiles


def get_hdf_datasets(fpath, pattern=None):
    """Read in an HDF4 (.hdf) file.

    Parameters
    ----------
    fpath : str
        File path to an HDF4 file.
    pattern : str | NoneType
        A pattern used to filter the names to only those containing it. This
        is case insensitive. Defaults to None, or no filtering.

    Returns
    -------
    list[str] : A list of strings representing all the variables names in the
        input HDF file.
    """
    # Get all variable identifier information in the file with GDAL
    ds = gdal.Open(fpath)
    datasets = ds.GetSubDatasets()

    # Convert to just a list of dataset names
    names = [d[0] for d in datasets]

    # If a pattern is given, filter these names for that
    if pattern:
        names = [n for n in names if pattern.lower() in n.lower()]

    return names


def hdf4_to_geotiff(fpath, dst, dataset=None, pattern=None):
    """Convert an HDF4 dataset to a GeoTiff.

    A dataset or pattern argument must be provided.

    Parameters
    ----------
    fpath : str
        File path to an HDF4 file.
    dataset : str
        The name of a variable in the HDF file.
    dst : str
        The target path for the output GeoTiff file.
    """
    # Open the HDF file and list all variable names
    ds = gdal.Open(fpath, gdal.GA_ReadOnly)
    names = [d[0] for d in ds.GetSubDatasets()]

    # Make sure a dataset or a pattern was given
    if not dataset and not pattern:
        raise ValueError(
            "No dataset name or dataset matching pattern provided for "
            f"{fpath}, cannot write file."
        )

    # Optionally filter these to match a pattern
    if not dataset and pattern:
        matches = [n for n in names if pattern.lower() in n.lower()]
        if len(matches) > 1:
            raise ValueError(
                "There are multiple datasets matching the pattern "
                f"'{pattern}' in {fpath}: {names}. Please choose a more "
                "precise pattern or provide the full name of a singular "
                f"dataset. Available datasets: {names}"
            )
        elif len(matches) == 0:
            raise ValueError(
                "There are no datasets matching the pattern "
                f"'{pattern}' in {fpath}. Please choose or precise "
                "pattern or provide the full name of a singular dataset. "
                "Available datasets: {names}"
            )
        dataset = matches[0]

    # Otherwise, make sure the target dataset is in this file
    else:
        if dataset not in names:
            raise ValueError(
                f"{dataset} not found in {fpath}, available datasets: {names}"
            )

    # We need to find the index position of the target variable
    idx = names.index(dataset)

    # Pull this one out into a Dataset object
    layer = gdal.Open(names[idx])

    # Warp this layer to the target file
    dst = os.path.expanduser(dst)
    gdal.Warp(dst, layer)
    del ds
