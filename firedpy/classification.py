# -*- coding: utf-8 -*-
"""Experimental methods for classifying burn data into discrete events.

Author: travis
Date: Sat Apr 25 12:38:15 PM MDT 2026
"""
import datetime as dt
import os

from functools import cached_property
from concurrent.futures import ProcessPoolExecutor, as_completed
from logging import getLogger
from pathlib import Path

import dask
import dask.array as da
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

from osgeo import gdal
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components
from shapely import Point
from sklearn.neighbors import radius_neighbors_graph
from tqdm import tqdm

gdal.UseExceptions()
logger = getLogger(__name__)


HOME = Path(__file__).parent


class Classifier:
    """Methods for classifying burn data into discrete fire events."""

    def __init__(self, hdf_dir, spatial_param=5, temporal_param=11,
                 shp_fpath=None, sample=False):
        """Initialize an Classifier object.

        Parameters
        ----------
        hdf_dir : str | pathlib.PosixPath
            Path to directory containing original MCD64A1 burn data HDF4 files
            with file extension '.hdf'. A nested file search will be performed
            on this directory and all files anywhere in it will be used.
        spatial_param : int
            The number of cells (~463 m resolution) to search for neighboring
            burn detections. Defaults to 5.
        temporal_param : int
            The number of days to search for neighboring burn detections.
            Defaults to 11.
        shp_fpath : str
            Path to a shapefile to use for the fire study area. Defaults to
            None.
        sample : bool
            Restrict classification to the first 10,000 files.
        """
        self.hdf_dir = hdf_dir
        self.spatial_param = spatial_param
        self.temporal_param = temporal_param
        self.shp_fpath = shp_fpath
        self.sample = sample

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
        r = gdal.Open(ds)
        geom = r.GetGeoTransform()
        df["x"] = df.apply(self.add_coordinate, geom=geom, which="x", axis=1)
        df["y"] = df.apply(self.add_coordinate, geom=geom, which="y", axis=1)
        return df

    def build_burns(self):
        """Build a composite data frame with just burn detections."""
        data = []
        with ProcessPoolExecutor(os.cpu_count()) as pool:
            jobs = [pool.submit(self.to_dataframe, fp) for fp in self.files]
            for job in tqdm(as_completed(jobs), total=len(jobs)):
                df = job.result()
                if df.shape[0] > 0:
                    data.append(df)

        # Create data frame
        df = pd.concat(data)
        df = df.sort_values("date")
        df = df.reset_index(drop=True)
        df["index"] = df.index

        return df

    @cached_property
    def shapefile(self):
        """Return a shapefile if a path is provided."""
        shp = None
        if self.shp_fpath:
            shp_name = Path(self.shp_fpath).name
            logger.info(f"Shapefile provided, events will be clipped to "
                        f"{shp_name}...")

            # Read in the shapefile
            shp = gpd.read_file(self.shp_fpath)
            shp.to_crs(self.profile["crs"], inplace=True)

            # Buffer to ensure fires immediately outside border are captured
            shp["geometry"] = shp["geometry"].buffer(100_000)
            shp.loc[:, "intersects"] = 1

        return shp

    def clip_to_shape(self, gdf):
        """Clip fire events to shapefile boundry, keep overlapping events.  

        Parameters
        ----------
        gdf : geopandas.geodataframe.GeoDataFrame
            GeoDataFrame of Firedpy fire events.
        shape_file : str | pathlib.PosixPath
            Path to shapefile representing a study area.
        buffer : int
            A distance in meters used to buffer the clipping shapefile to
            ensure fires immediately outside the border are captured.

        Returns
        -------
        geopandas.geodataframe.GeoDataFrame : A clipped GeoDataFrame.
        """
        # Characterize shapefile intersecting events by ID  # <---------------- Clipping earlier at parallel read step, but that precludes classification which is required here (check that this is acceptable)
        gdf = gpd.sjoin(gdf, self.shapefile, how="left")
        gdf.loc[:, "keep"] = gdf.groupby("id")["intersects"].transform("any")

        # Drop non-intersecting events
        clipped_gdf = gdf[gdf["keep"]]
        for tmp_field in ["index_right", "intersects", "keep"]:
            del clipped_gdf[tmp_field]

        return clipped_gdf

    @property
    def files(self):
        """Return all HDF4 files in given directory."""
        files = list(Path(self.hdf_dir).rglob("*hdf"))
        if self.sample:
            files = files[:1_000]
        return files

    @property
    def profile(self):
        """Return geographic attributes of an HDF4 file."""
        sample_fpath = self.files[0]
        ds = gdal.Open(sample_fpath).GetSubDatasets()[0][0]
        r = gdal.Open(ds)
        x = r.RasterXSize
        y = r.RasterYSize
        gtype = r.GetRasterBand(1).DataType
        geom = r.GetGeoTransform()
        attrs = {
            "dtype":gdal.GetDataTypeName(gtype).lower(),
            "shape": (y, x),
            "transform": geom,
            "resolution": geom[1],
            "crs": r.GetProjection()
        }
        return attrs

    @dask.delayed
    def read_band(self, h4_fpath):
        """Read the HDF4 Burn Date array as a Dask Delayed object."""
        ds = gdal.Open(h4_fpath).GetSubDatasets()[0][0]
        r = gdal.Open(ds)
        obj = r.GetRasterBand(1).ReadAsArray()
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
        dtype = self.profile["dtype"]
        shape = self.profile["shape"]

        # Convert to a dask data frame
        data = da.from_delayed(ddata, shape=shape, dtype=dtype)
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
        df = self.to_geo(df)

        return df

    def to_geo(self, df):
        """Convert dataframe to geodataframe and clip if shapefile provided.

        Parameters
        ----------
        df : pandas.core.frame.DataFrame
            A pandas dataframe with x, y coordinates.

        Returns
        -------
        geopandas.geodataframe.GeoDataFrame
            A GeoDataFrame of the original DataFrame object.
        """
        # Build a datetime coordinate dataframe from the events
        logger.info("Converting burn detection dataframe to geodataframe...")

        # Center pixel coordinates
        df.loc[:, "x"] = df["x"] + (self.profile["transform"][1] / 2)
        df.loc[:, "y"] = df["y"] + (self.profile["transform"][-1] / 2)

        # Each entry gets a point object from the x and y coordinates.
        df.loc[:, "geometry"] = df[["x", "y"]].apply(
            lambda x: Point(tuple(x)),
            axis=1
        )
        df = gpd.GeoDataFrame(df, crs=self.profile["crs"], geometry="geometry")

        if self.shp_fpath:
            df = gpd.clip(df, self.shapefile)

        return df

    def classify(self):
        """Classify burn detections into discrete events."""
        # Get the burn detection dataframe
        logger.info("Converting burn detections to data frame...")
        df = self.build_burns()

        # If we scale the space and time coordinates, we can use a radius of 1
        coords = df[["x", "y", "day"]].values.astype("float32")
        coords[:, 0:2] /= self.spatial_param * self.profile["resolution"]
        coords[:, 2] /= self.temporal_param

        # This creates a compressed sparse graph from the coordinates
        graph = radius_neighbors_graph(
            coords, 
            radius=1.0, 
            metric="chebyshev",  # Also called "chessboard distance"
            mode="connectivity", 
            n_jobs=-1
        )

        # Connects each observation if they are within 1 in all directions
        n_components, labels = connected_components(
            csgraph=graph, 
            directed=False, 
            return_labels=True
        )
        df["id"] = labels

        return df


if __name__ == "__main__":
    self = Classifier(
        hdf_dir=Path("/home/travis/scratch/firedpy/conus_2000_2026"),
        sample=True
    )
    df = self.classify()
    print(df)