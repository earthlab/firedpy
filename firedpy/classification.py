# -*- coding: utf-8 -*-
"""Experimental methods for classifying burn data into discrete events.

Author: travis
Date: Sat Apr 25 12:38:15 PM MDT 2026
"""
import datetime as dt

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import dask
import dask.array as da
import numpy as np
import pandas as pd
import xarray as xr

from osgeo import gdal
from tqdm import tqdm

gdal.UseExceptions()


HOME = Path(__file__).parent


class Classifier:
    """Methods for classifying burn data into discrete fire events."""

    def __init__(self, hdf_dir, spatial_param=5, temporal_param=11):
        """Initialize an Classifier object."""
        self.hdf_dir = hdf_dir
        self.spatial_param = spatial_param
        self.temporal_param = temporal_param

    def __repr__(self):
        name = self.__class__.__name__
        args = {key: self.__dict__[key] for key in self.__static_attributes__}
        address = hex(id(self))
        msgs = [f"\n   {k}={v}" for k, v in args.items()]
        msg = " ".join(msgs)
        return f"<{name} object at {address}> {msg}"

    def add_coordinate(self, row, geom, which):
        """Convert index positions in a pandas series to coordinates."""
        if which == "x":
            c = row["x"] * geom[1] + geom[0]
        elif which == "y":
            c = row["y"] * geom[5] + geom[3]
        else:
            raise ValueError("`which` parameter must be 'x' or 'y', received "
                             f"'{which}.")
        return c

    def add_coordinates(self, df, h4_fpath):
        """Add coordinates using x, y positions and the original file."""
        ds = gdal.Open(h4_fpath).GetSubDatasets()[0][0]
        ds = gdal.Open(ds)
        geom = ds.GetGeoTransform()
        df["x"] = df.apply(self.add_coordinate, geom=geom, which="x", axis=1)
        df["y"] = df.apply(self.add_coordinate, geom=geom, which="y", axis=1)
        return df

    @property
    def files(self):
        """Return all HDF4 files in given directory."""
        files = list(Path(self.hdf_dir).rglob("*hdf"))
        return files

    def get_attrs(self, h4_fpath):
        """Return specific attributes of an HDF4 file."""
        ds = gdal.Open(h4_fpath).GetSubDatasets()[0][0]
        ds = gdal.Open(ds)
        x = ds.RasterXSize
        y = ds.RasterYSize
        gtype = ds.GetRasterBand(1).DataType
        attrs = {
            "dtype":gdal.GetDataTypeName(gtype).lower(),
            "shape": (y, x)
        }
        return attrs

    @dask.delayed
    def read_band(self, h4_fpath):
        """Read the HDF4 Burn Date array as a Dask Delayed object."""
        ds = gdal.Open(h4_fpath).GetSubDatasets()[0][0]
        ds = gdal.Open(ds)
        obj = ds.GetRasterBand(1).ReadAsArray()
        return obj

    def to_date(self, doy, h4_fpath, since=1970):
        """Convert a day of year value to days since a given year."""
        date_string = h4_fpath.name.split(".")[1]
        year = int(date_string[1:5])
        date = dt.datetime(year, 1, 1) + dt.timedelta(days=int(doy) - 1)
        return date

    def to_dataframe(self, h4_fpath=None):
        """Convert the original MODIS HDF4 dataset into a data frame.

        NOTE: There are pithier Xarray methods that can read this data in
            lazily, but the backend engines appear to require non-pypi binaires
            or else aren't ready for Python 3.14.
        """
        # Get delayed data, shape, and data type attributes of file
        ddata = self.read_band(h4_fpath)
        attrs = self.get_attrs(h4_fpath)

        # Convert to a dask data frame
        data = da.from_delayed(ddata, **attrs)
        darray = xr.DataArray(data, dims=("y", "x"), name="value")
        ddf = darray.to_dask_dataframe()

        # Filter for valid dates and pull into memory
        ddf = ddf[ddf["value"] > 0]
        df = ddf.compute()

        # Convert date
        df["date"] = df["value"].apply(self.to_date, h4_fpath=h4_fpath)

        # Conver this to days since 1970 for a nice integer
        base = dt.datetime(1970, 1, 1)
        df["day"] = df["date"].apply(lambda d: (d - base).days)

        # Add the tile ID, just in case it's useful
        df["tile"] = h4_fpath.name.split(".")[2]

        # The x, y coordinates here are just the array's index position
        df = self.add_coordinates(df, h4_fpath)

        return df

    def build_burns(self):
        """Build a composite data frame with just burn detections."""
        data = []
        with ThreadPoolExecutor() as pool:
            jobs = [pool.submit(self.to_dataframe, fp) for fp in self.files]
            for job in tqdm(jobs):
                df = job.result()
                if df.shape[0] > 0:
                    data.append(df)
        df = pd.concat(data)
        df = df.sort_values("date")
        df = df.reset_index(drop=True)
        df["index"] = df.index
        return df

    def classify(self):
        """Classify burn detections into discrete events."""
        # Experimenting with a "connected components" algorithm
        from scipy.sparse.csgraph import connected_components
        from scipy.sparse import csr_array
        from sklearn.neighbors import radius_neighbors_graph

        # Get the burn detection dataframe
        df = self.build_burns()  # For the US from 2014-2026, this takes about 2 and half minutes and maxes out at about 12.75 GB on my machine (including 5 GB of overhead)

        # If we scale the space and time coordinates, we can use a radius of 1
        coords = df[["x", "y", "day"]].values.astype("float32")
        coords[:, 0:2] /= self.spatial_param * 500  # 500m is a standin for the pixel size
        coords[:, 2] /= self.temporal_param

        # This creates a compressed sparse graph from the coordinates
        graph = radius_neighbors_graph(
            coords, 
            radius=1.0, 
            metric="chebyshev",  # Also called "chessboard distance"
            mode="connectivity", 
            n_jobs=-1
        )

        # And this connects each observation if they are within 1 in all
        # directions
        n_components, labels = connected_components(
            csgraph=graph, 
            directed=False, 
            return_labels=True
        )
        df["event"] = labels

        return df

    def to_geo(self, tdf):
        """Convert dataframe to geodataframe."""
        from shapely import Point
        import geopandas as gpd

        tdf["geometry"] = tdf.apply(lambda p: Point(p["x"], p["y"]), axis=1)
        proj = 'PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid",DATUM["Not specified (based on custom spheroid)",SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'
        tdf = gpd.GeoDataFrame(tdf, geometry="geometry", crs=proj)
        tdf.to_parquet("~/scratch/fried_test.parquet")


if __name__ == "__main__":
    hdf_dir = Path("/home/travis/scratch/firedpy/conus_2014_2026/rasters/burn_area/hdfs/")
    self = Classifier(hdf_dir)
