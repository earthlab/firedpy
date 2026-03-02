import hashlib
import math
import random
import re

from argparse import Namespace
from concurrent.futures import as_completed, ProcessPoolExecutor
from datetime import datetime
from logging import getLogger
from pathlib import Path, PosixPath
from typing import List, Set, Tuple

import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
import xarray as xr

from pyproj import CRS, Proj, transform
from rasterio.features import rasterize
from rasterio.transform import from_bounds
from scipy.spatial import cKDTree
from shapely import Point, Polygon, MultiPolygon
from tqdm import tqdm

from firedpy import DATA_DIR
from firedpy.data_classes import Base
from firedpy.utilities.spatial import MODIS_CRS, get_country_file

logger = getLogger(__name__)


def generate_path(project_directory, base_filename, shape_type):
    """Return the full paths to a target files.

    Parameters
    ----------
    project_directory : str
        A project directory path containing target shapefiles. This
        corresponds with the '-o' or '--output_directory` option in the
        firedpy CLI.
    base_filename : str
        The file name of the target file.
    shape_type : str
        Build shapefiles from the event data frame. Specify either "shp",
        "gpkg", or both. Shapefiles of both daily progression and overall
        event perimeters will be written to the 'outputs/shapefiles' folder of
        the chosen project directory. These will be saved in the specified
        geopackage format (.gpkg), ERSI Shapefile format (.shp), or save them
        in both formats using the file basename of the fire event data frame
        (e.g. 'modis_events_daily.gpkg' and 'modis_events.gpkg').

    Returns
    -------
    list[str] : A list of full filepaths corresponding with the target file
        and file formats.
    """
    # Get the requested combination of shapefile types
    file_extensions = {
        "shp": [".shp", None],
        "gpkg": [None, ".gpkg"],
        "both": [".shp", ".gpkg"]
    }
    file_ext = file_extensions.get(shape_type, [".gpkg"])

    # Create the paths
    paths = []
    proj_dir = Path(project_directory)
    proj_dir = proj_dir.expanduser().absolute()
    for ext in file_ext:
        if ext:
            fname = f"{base_filename}{ext}"
            path = proj_dir.joinpath("outputs", "shapefiles", fname)
        else:
            path = None
        paths.append(path)

    return paths


def _merge_events_spatially_task(namespace: Namespace):
    namespace, pos = namespace
    events: List[EventPerimeter] = namespace.events
    sp_buf: float = namespace.sp_buf

    if len(events) == 1:
        return events

    ckd_trees = []
    for event in events:
        coords = np.array([(c[1], c[0]) for c in event.spacetime_coordinates])
        ckd_trees.append(cKDTree(coords))

    i, j = 0, 1
    merged_event_indices = set()
    merged_events = []
    pbar = tqdm(
        total=len(ckd_trees),
        desc=f"Merging edges for temporal group {pos}",
        position=pos
    )
    found_overlapping = False
    while i < len(ckd_trees):
        if j not in merged_event_indices and j < len(ckd_trees):
            indices = ckd_trees[i].query_ball_tree(ckd_trees[j], r=sp_buf, p=2)
            if any(index_list for index_list in indices):
                found_overlapping = True
                coords = events[j].spacetime_coordinates
                events[i].add_spacetime_coordinates(coords)
                events[i].is_edge = events[i].is_edge or events[j].is_edge
                crds = [(c[1], c[0]) for c in events[i].spacetime_coordinates]
                ckd_trees[i] = cKDTree(np.array(crds))
                merged_event_indices.add(j)

        j += 1

        if j >= len(ckd_trees):
            if not found_overlapping:
                merged_events.append(events[i])
                i += 1
                while i in merged_event_indices:
                    i += 1

            pbar.update(i)
            j = i + 1

            found_overlapping = False

    return merged_events


class SpacetimeCoordinate:
    __slots__ = ('x', 'y', 't')  # Helps to save memory

    def __init__(self, x: float, y: float, t: int):
        # Ensure t within 16-bit unsigned integer range
        if not (0 <= t <= 65535):
            raise ValueError("t must be a 16-bit unsigned integer (0-65535)")

        self.x = np.float64(x)
        self.y = np.float64(y)
        self.t = np.uint16(t)

    def __repr__(self):
        return f"SpacetimeCoordinate(x={self.x}, y={self.y}, t={self.t})"

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.t == other.t

    def __hash__(self):
        # To make object hashable, return hash based on values of x, y, and t
        return hash((self.x, self.y, self.t))

    # The following methods are added for pickling support
    def __getstate__(self):
        return self.x, self.y, self.t

    def __setstate__(self, state):
        self.x, self.y, self.t = state


class EventPerimeter:
    def __init__(
            self,
            event_id: int,
            spacetime_coordinates: Set[Tuple[float, float, np.int16]],
            is_edge: bool = False
    ):
        self.event_id = event_id
        self.spacetime_coordinates = spacetime_coordinates
        self.is_edge = is_edge
        self.min_y = self.max_y = self.min_x = self.max_x = None
        self.min_t = self.max_t = None
        self.min_geom_x = None
        self.max_geom_y = None

    def add_spacetime_coordinates(
            self,
            new_coordinates: Set[Tuple[float, float, np.int16]]
    ):
        self.spacetime_coordinates.update(new_coordinates)

    def compute_min_max(self):
        # Convert list of tuples to a NumPy array for efficient processing
        array = np.array(list(self.spacetime_coordinates))

        # Check that the array is of the correct shape
        if array.ndim == 2 and array.shape[1] == 3:
            # Find min and max for each column (dimension)
            self.min_y, self.min_x, self.min_t = array.min(axis=0)
            self.max_y, self.max_x, self.max_t = array.max(axis=0)

            # Return the min and max values as a convenience
            return {
                "min_y": self.min_y,
                "min_x": self.min_x,
                "min_t": self.min_t,
                "max_y": self.max_y,
                "max_x": self.max_x,
                "max_t": self.max_t
            }
        else:
            raise ValueError("Array should be of shape (N, 3) where N is the "
                             "number of 3-tuples.")

    def __add__(self, other):
        if not isinstance(other, EventPerimeter):
            raise TypeError('Can only add other EventPerimeter')
        self.spacetime_coordinates.update(other.spacetime_coordinates)
        return self

    def __eq__(self, other):
        if not isinstance(other, EventPerimeter):
            return NotImplemented
        return self.event_id == other.event_id

    def __hash__(self):
        return hash(self.event_id)


class RandomAccessSet:
    def __init__(self):
        self.list = []
        self.set = {}

    def add(self, value):
        if value not in self.set:
            self.list.append(value)
            self.set[value] = len(self.list) - 1

    def remove(self, value):
        # Remove in constant time
        if value in self.set:
            index_to_remove = self.set[value]
            last_element = self.list[-1]
            self.list[index_to_remove], self.list[-1] = last_element, value
            self.set[last_element] = index_to_remove
            self.list.pop()
            del self.set[value]

    def get_random(self):
        if not self.list:
            raise ValueError("RandomAccessSet is empty")
        return random.choice(self.list)

    def get(self, value):
        if value in self.set:
            return self.list[self.set[value]]


class EventGrid(Base):
    """Gridded event delineation methods.

    For a single file, however large, find sites with any burn detections
    in the study period, loop through these sites and group by the space-time
    window, save grouped events to a data frame on disk.
    """

    def __init__(
            self,
            project_directory,
            spatial_param=5,
            temporal_param=11,
            start_year=2020,
            end_year=2025,
            country=None,
            shape_file=None,
            area_unit="Unknown",
            time_unit="days since 1970-01-01",
            nc_fpath=None,
            input_array=None,
            coordinates=None
    ):
        """Initialize EventGrid object.

        TODO: Check the definitions/types of the last two arguments here.

        Parameter
        ---------
        project_directory : str | pathlib.PosixPath
            Path to firedpy project directory.
        spatial_param : int
            The number of cells (~463 m resolution) to search for neighboring
            burn detections. Defaults to 5 cells in all directions. Defaults
            to 5.
        temporal_param : int
            The number of days to search for neighboring burn detections.
            Defaults to 11.
        start_year : int
            The first year of fire events. Defaults to 2000.
        end_year : int
            The last year of fire events. Defaults to 2025.
        country : str
            The name of a country to use as a study area.
        shape_file : str
            Path to a shapefile to use for the fire study area. Defaults to
            None.
        area_unit : str
            Area Unit. Currently unused and defaults to 'Unknown'.
        time_unit : str
            The time unit to use when parsing dates. Defaults to days since
            1970-01-01.
        nc_fpath : str | pathlib.PosixPath
            Path to the input NetCDF4 file containing the original burn data.
        input_array : np.ndarray
        coordinates : dict[str, np.ndarray]
        """
        super().__init__(project_directory)
        self.spatial_param = spatial_param
        self.temporal_param = temporal_param
        self.shape_file = shape_file
        self.current_event_id = 0
        self.years = list(range(start_year, end_year + 1))
        if country:
            self.shape_file = get_country_file(country)
        elif shape_file:
            self.shape_file = shape_file

        self._area_unit = area_unit
        self._time_unit = time_unit

        # If a country was provided, convert to a shapefile
        if nc_fpath:
            self.input_array = self.build_array(nc_fpath)

        elif input_array is not None:
            if coordinates is None:
                msg = "Coordinates parameter also required with input array."
                raise ValueError(msg)
            self._coordinates = coordinates
            self.input_array = input_array

        else:
            msg = "Must input either nc_file_path or input_array parameter"
            raise ValueError(msg)

        del self.n_cores

    def __repr__(self):
        """Return representation string for BurnData object."""
        name = self.__class__.__name__
        address = hex(id(self))
        attrs = {}
        msgs = []
        for key, attr in self.__dict__.items():
            if "data_frame" in key:
                attr = f"{type(attr)} {attr.shape}"  # Too big for preview
            if not key.startswith("_"):  # Avoid secrets/private attributes
                attrs[key] = attr
            if key == "input_array":
                attrs[key] = f"{type(attr)}: {attr.shape}"
        for key, value in attrs.items():
            if isinstance(value, (str, PosixPath)):
                msgs.append(f"\n   {key}='{value}'")
            else:
                msgs.append(f"\n   {key}={value}")
        msg = " ".join(msgs)
        return f"<{name} object at {address}> {msg}"

    def build_array(self, nc_fpath):
        """Read and subset the array to be converted to fire events values.

        Parameters
        ---------
        nc_fpath : str | pathlib.PosixPath
            A path to a NetCDF file of MODIS burn dates with "value" as the
            date dataset.

        Returns
        -------
        numpy.ndarray : A numpy array.
        """
        # Read in the netcdf data array
        with xr.open_dataset(nc_fpath) as burns:
            burns = burns.sel(time=burns.time.dt.year.isin(self.years))
            self._coordinates = {}
            if self.shape_file:
                burns = self.mask_array(burns, self.shape_file)
            for dim in ["y", "x", "time"]:
                coords = np.array(burns.coords[dim].values)
                self._coordinates[dim] = coords
        return burns.value.values

    def mask_array(self, array, shape_file):
        """Mask an xarray dataset with a shapefile.

        Parameters
        ----------
        array : xarray.core.dataset.Dataset
            An xarray dataset.
        shape_file : str | pathlib.PosixPath
            Path to a shapefile to use to mask the array with nans.

        Returns
        -------
        xarray.core.dataset.Dataset : A masked xarray dataset.
        """
        # Check that the coordinate reference systems match
        mask = gpd.read_file(shape_file)
        target_crs = CRS(array.crs.spatial_ref).to_wkt()
        if mask.crs != target_crs:
            mask = mask.to_crs(target_crs)

        # Buffer mask to ensure fires immediately outside border are captured
        mask["geometry"] = mask["geometry"].buffer(100_000)

        # Get the geometries used to burn values to the array
        shapes = ((geom, 1) for geom in mask["geometry"])

        # Get the target array geometry
        ny, nx = array.sizes["y"], array.sizes["x"]
        xmin, xmax = array.x.values[0], array.x.values[-1]
        ymax, ymin = array.y.values[0], array.y.values[-1]
        transform = from_bounds(xmin, ymin, xmax, ymax, nx, ny)

        # Rasterize mask
        mask_array = rasterize(
            shapes=shapes,
            out_shape=(ny, nx),
            transform=transform,
            fill=0,
            all_touched=False
        )

        # mask original array
        array["value"] *= mask_array

        return array

    def _get_spatial_window(self, y, x, array_dims):
        """
        Pull in the spatial window around a detected event and determine its
        shape and the position of the original point within it. Finding this
        origin point is related to another time-saving step in the event
        classification procedure.
        """
        top = max(0, y - self.spatial_param)
        bottom = min(array_dims[0], y + self.spatial_param)
        left = max(0, x - self.spatial_param)
        right = min(array_dims[1], x + self.spatial_param)

        def get_edge_positions(dim, coord, param):
            """Compute center and origin."""
            edges_lower = [t for t in range(param)]
            edges_upper = [dim - t - 1 for t in range(param)]
            is_edge = True
            if coord in edges_lower:
                center = coord
                origin = 0
            elif coord in edges_upper:
                center = coord - (coord - param if coord >= param else 0)
                origin = coord - param if coord >= param else 0
            else:
                is_edge = False
                center = param
                origin = coord - param  # top left

            return center, origin, is_edge

        dim0 = array_dims[0]
        dim1 = array_dims[1]
        sp = self.spatial_param
        ycenter, oy, is_y_edge = get_edge_positions(dim0, y, sp)
        xcenter, ox, is_x_edge = get_edge_positions(dim1, x, sp)

        edge = is_x_edge or is_y_edge

        return top, bottom, left, right, [ycenter, xcenter], [oy, ox], edge

    def _get_available_cells(self):
        """
        To save time, avoid checking cells with no events at any time step.
        Create a mask of max values at each point. If the maximum at a cell is
        less than or equal to zero there were no values and it will not be
        checked in the event classification step.
        """
        # Using memory - can handle large tiles, but gets pretty high
        locs = np.where(np.any(self.input_array > 0, axis=0))

        return list(zip(locs[0], locs[1]))

    def get_event_perimeters(self, all_t=False):
        """Group MODIS burn detections cells into fire events.

        Parameters
        ----------
        all_t : bool
            All touch.

        Returns
        -------
        list[firedpy.model_classes.EventPerimeter] : A list of fire event
            perimeter objects.
        """
        available_pairs = self._get_available_cells()
        nz, ny, nx = self.input_array.shape
        dims = [ny, nx]
        perimeters = RandomAccessSet()
        identified_points = {}

        time_index_buffer = max(1, self.temporal_param // 30)

        for pair in available_pairs:
            y, x = pair

            top, bottom, left, right, center, origin, is_edge = \
                self._get_spatial_window(y, x, dims)

            center_y, center_x = center
            window = self.input_array[:, top:bottom + 1, left:right + 1]
            center_burn = window[:, center_y, center_x]
            valid_center_burn_indices = np.where(center_burn > 0)[0]
            valid_center_burn_values = center_burn[valid_center_burn_indices]

            for i, burn_value in enumerate(valid_center_burn_values):
                overlapping_event_ids = RandomAccessSet()
                events_to_remove = set()

                # In each time slice we have a possible range of 30 days
                if all_t:
                    left = 0
                    right = nz
                else:
                    n = int(valid_center_burn_indices[i] - time_index_buffer)
                    left = max(0, n)
                    right = min(nz, n + 1)

                burn_window = window[left:right, :, :]
                val_locs = np.where(
                    abs(burn_value - burn_window) <= self.temporal_param
                )
                y_locs, x_locs = val_locs[1], val_locs[2]
                oy, ox = origin

                vals = burn_window[val_locs]

                ys = self._coordinates["y"][oy + y_locs]
                xs = self._coordinates["x"][ox + x_locs]

                new_spacetime_coords = set()
                for i, val in enumerate(vals):
                    c = (np.float32(ys[i]), np.float32(xs[i]), np.uint16(val))

                    if c in identified_points:
                        overlapping_event_ids.add(identified_points[c])
                    else:
                        new_spacetime_coords.add(c)

                event = EventPerimeter(
                    event_id=self.current_event_id,
                    spacetime_coordinates=new_spacetime_coords,
                    is_edge=is_edge
                )

                if overlapping_event_ids.list:
                    event_ids = [e for e in overlapping_event_ids.list]
                    to_keep_ids = sorted(
                        event_ids,
                        key=lambda event: len(
                            perimeters.get(
                                EventPerimeter(event, set())
                            ).spacetime_coordinates
                        )
                    )
                    to_keep_id = to_keep_ids[-1]
                    to_keep = perimeters.get(EventPerimeter(to_keep_id, set()))
                    for event_id in overlapping_event_ids.set:
                        if event_id != to_keep_id:
                            to_merge = perimeters.get(
                                EventPerimeter(event_id, set())
                            )
                            to_keep.add_spacetime_coordinates(
                                to_merge.spacetime_coordinates
                            )
                            is_edge = to_keep.is_edge or to_merge.is_edge
                            to_keep.is_edge = is_edge
                            coords = to_merge.spacetime_coordinates
                            identified_points.update(
                                (p, to_keep_id) for p in coords
                            )
                            events_to_remove.add(to_merge)

                    to_keep.add_spacetime_coordinates(
                        event.spacetime_coordinates
                    )
                    to_keep.is_edge = to_keep.is_edge or event.is_edge
                    identified_points.update(
                        (p, to_keep_id) for p in new_spacetime_coords
                    )

                    if to_keep.is_edge:
                        if to_keep.min_geom_x is None:
                            to_keep.min_geom_x = self._coordinates["x"][0]
                        if to_keep.max_geom_y is None:
                            to_keep.max_geom_y = self._coordinates["y"][0]

                elif event.spacetime_coordinates:
                    if event.is_edge:
                        event.min_geom_x = self._coordinates["x"][0]
                        event.max_geom_y = self._coordinates["y"][0]

                    perimeters.add(event)
                    identified_points.update(
                        (p, event.event_id) for p in new_spacetime_coords
                    )
                    self.current_event_id += 1

                for event in events_to_remove:
                    perimeters.remove(event)

        return perimeters.list


class ModelBuilder(Base):
    def __init__(
            self,
            project_directory,
            tiles,
            country=None,
            shape_file=None,
            spatial_param=5,
            temporal_param=11,
            start_year=2020,
            end_year=2025,
            n_cores=0
    ):
        """Methods for classifying fire events.

        Parameters
        ----------
        project_directory : str
            Project output directory path. Required.
        tiles : str | list
            A string representing a single MODIS tile (e.g., 'h08v04'), a
            string representing multiple tiles separated by spaces
            (e.g., 'h08v04 h09v04') or a list representing multiple tiles
            (e.g., ['h08v04', 'h09v04']).
        country : str
            The name of a country to use as a study area. If not provided,
            either 'tiles' or 'shape_file' parameter must be provided.
        shape_file : str
            Path to a shapefile to use for the fire study area. Defaults to
            None.
        spatial_param : int
            The number of cells (~463 m resolution) to search for neighboring
            burn detections.
        temporal_param : int
            The number of days to search for neighboring burn detections.
        start_year : int
            The first year of fire events. Defaults to 2000.
        end_year : int
            The last year of fire events. Defaults to 2025.
        n_cores : int
            Number of cores to use for parallel processing. A value of 0 or
            None will use all available cores. Defaults to 0.
        """
        super().__init__(project_directory, n_cores)
        self.tiles = tiles
        self.country = country
        self.shape_file = shape_file
        self.spatial_param = spatial_param
        self.temporal_param = temporal_param
        self.start_year = start_year
        self.end_year = end_year
        self._lc_mosaic_re = r'lc_mosaic_(?P<land_cover_type>\d{1})_\
            (?P<year>\d{4})\.tif$'

        # Use the first file to get some geometry data for later
        with xr.open_dataset(self.files[0]) as data_set:
            self.crs = data_set.crs
            self.geom = self.crs.geo_transform
            self._res = self.geom[1]
            self.sp_buf = spatial_param * self._res
            dims = ["y", "x", "time"]
            self._coordinates = {
                dim: np.array(data_set.coords[dim].values) for dim in dims
            }

    def __repr__(self):
        """Return representation string for BurnData object."""
        name = self.__class__.__name__
        address = hex(id(self))
        attrs = {}
        msgs = []
        for key, attr in self.__dict__.items():
            if "data_frame" in key:
                attr = f"{type(attr)} {attr.shape}"  # Too big for preview
            if not key.startswith("_"):  # Avoid secrets/private attributes
                attrs[key] = attr
            if key == "crs":
                attrs[key] = CRS(attr.spatial_ref).to_proj4()
        for key, value in attrs.items():
            if isinstance(value, (str, PosixPath)):
                msgs.append(f"\n   {key}='{value}'")
            else:
                msgs.append(f"\n   {key}={value}")
        msg = " ".join(msgs)
        return f"<{name} object at {address}> {msg}"

    def add_fire_attributes(self, gdf):
        """Add fired attributes to an event geodataframe.

        Parameters
        ----------
        gdf : geopandas.geodataframe.GeoDataFrame
            A geodataframe of fire events.

        Returns
        -------
        geopandas.geodataframe.GeoDataFrame : The geodataframe with fire
            attribute fields appended.
        """
        if gdf.shape[0] == 0:
            logger.warning("No fire events found")
            return gdf

        logger.info("Adding fire attributes ...")
        gdf['pixels'] = gdf.groupby(['id', 'date'])['id'].transform('count')
        gdf['ig_utm_x'] = gdf.groupby(['id', 'date'])['x'].nth(0)
        gdf['ig_utm_y'] = gdf.groupby(['id', 'date'])['y'].nth(0)
        group = gdf.groupby('id')
        gdf['date'] = gdf['date'].apply(
            lambda x: datetime.strptime(x, '%Y-%m-%d')
        )

        gdf['ig_date'] = group['date'].transform('min')
        gdf['ig_day'] = gdf['ig_date'].apply(
            lambda x: datetime.strftime(x, '%j')
        )

        gdf['ig_month'] = gdf['ig_date'].apply(lambda x: x.month)
        gdf['ig_year'] = gdf['ig_date'].apply(lambda x: x.year)
        gdf['last_date'] = group['date'].transform('max')

        gdf['tot_pix'] = group['id'].transform('count')

        gdf['daily_duration'] = gdf['date'] - gdf['ig_date']
        gdf['event_day'] = gdf['daily_duration'].apply(lambda x: x.days + 1)
        gdf['final_duration'] = gdf['last_date'] - gdf['ig_date']
        gdf['event_dur'] = gdf['final_duration'].apply(lambda x: x.days + 1)
        gdf.drop('final_duration', axis=1)

        gdf['dy_ar_km2'] = gdf['pixels'].apply(self._to_kms)
        gdf['tot_ar_km2'] = gdf['tot_pix'].apply(self._to_kms)

        gdf['fsr_px_dy'] = gdf['tot_pix'] / gdf['event_dur']
        gdf['fsr_km2_dy'] = gdf['fsr_px_dy'].apply(self._to_kms)

        gdf['mx_grw_px'] = group['pixels'].transform('max')
        gdf['mn_grw_px'] = group['pixels'].transform('min')
        gdf['mu_grw_px'] = group['pixels'].transform('mean')

        gdf['mx_grw_km2'] = gdf['mx_grw_px'].apply(self._to_kms)
        gdf['mn_grw_km2'] = gdf['mn_grw_px'].apply(self._to_kms)
        gdf['mu_grw_km2'] = gdf['mu_grw_px'].apply(self._to_kms)

        max_date = pd.DataFrame(group[['date', 'pixels']].apply(
            self._max_growth_date).reset_index()
        )
        max_date = max_date.rename(columns={0: 'mx_grw_dte'})
        gdf = gdf.merge(max_date[['id', 'mx_grw_dte']], on="id")

        gdf = gdf[['id', 'date', 'ig_date', 'ig_day', 'ig_month',
                   'ig_year', 'last_date', 'event_day', 'event_dur',
                   'pixels', 'tot_pix', 'dy_ar_km2', 'tot_ar_km2',
                   'fsr_px_dy', 'fsr_km2_dy', 'mx_grw_px', 'mn_grw_px',
                   'mu_grw_px', 'mx_grw_km2', 'mn_grw_km2', 'mu_grw_km2',
                   'mx_grw_dte', 'x', 'y', 'geometry', 'ig_utm_x',
                   'ig_utm_y']]

        gdf = gdf.reset_index(drop=True)

        return gdf

    ## Note from Nate: this is the definition that was exporting the data as WGS84
    ## Need to re-work this attribute type if we want to keep it
    ## currently commented out where called in build_events (line 855)
    def add_kg_attributes(self, gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """Assign Köppen-Geiger climate zones to events with a raster file.

        Args:
            gdf: GeoDataFrame containing 'x', 'y', and 'id'.
            tif_path: File path to a Köppen-Geiger .tif raster.

        Returns:
            GeoDataFrame with kg_zone, kg_mode, and kg_desc columns added.
        """
        tif_path = DATA_DIR.joinpath(
            'koppen_geiger_tif', '1991_2020', 'koppen_geiger_0p00833333.tif'
        )

        KG_LEGEND = {
            1: "Af", 2: "Am", 3: "Aw", 4: "BWh", 5: "BWk", 6: "BSh", 7: "BSk",
            8: "Csa", 9: "Csb", 10: "Csc", 11: "Cwa", 12: "Cwb", 13: "Cwc",
            14: "Cfa", 15: "Cfb", 16: "Cfc", 17: "Dsa", 18: "Dsb", 19: "Dsc",
            20: "Dsd", 21: "Dwa", 22: "Dwb", 23: "Dwc", 24: "Dwd", 25: "Dfa",
            26: "Dfb", 27: "Dfc", 28: "Dfd", 29: "ET", 30: "EF",
        }

        sgdf = gdf.copy()

        with rasterio.open(tif_path) as src:
            # Project GeoDataFrame to match raster CRS
            if sgdf.crs != src.crs:
                sgdf = sgdf.to_crs(src.crs)

            # Extract centroid coordinates from geometry (handles points,
            # polygons, multipolygons)
            coords = [
                (geom.centroid.x, geom.centroid.y) for geom in sgdf.geometry
            ]

            sampled_vals = list(src.sample(coords))
            sgdf['kg_zone'] = [
                v[0] if v[0] != src.nodata else np.nan for v in sampled_vals
            ]

        sgdf['kg_mode'] = sgdf.groupby('id')['kg_zone'].transform(self._mode)
        sgdf['kg_desc'] = sgdf['kg_mode'].map(KG_LEGEND)

        return sgdf

    def adjust_for_esri(self, gdf):
        """Adjust a geodataframe for the ESRI Shapefile driver."""
        sdf = gdf.copy()
        if "date" in sdf.columns:
            sdf["date"] = [str(d) for d in sdf["date"]]
        sdf["ig_date"] = [str(d) for d in sdf["ig_date"]]
        sdf["last_date"] = [str(d) for d in sdf["last_date"]]
        sdf["mx_grw_dte"] = [str(d) for d in sdf["mx_grw_dte"]]
        return sdf

    @staticmethod
    def _as_multi_polygon(polygon):
        if isinstance(polygon, Polygon):
            polygon = MultiPolygon([polygon])
        return polygon

    def build_events(self):
        """Build the fire event geodataframe."""
        # Classify MODIS burn data into fire events
        event_perimeters = self.classify_events()

        # Build the event point geodataframe
        gdf = self.build_points(event_perimeters)
 
        # Add fire event attributes
        gdf = self.add_fire_attributes(gdf)

        # Convert points to pixels
        gdf = self.process_geometry(gdf)

        # Calculate fire spread speed and maximum travel vectors
        # Note from Nate: this does not calculate fire speed and max travel vectors? And we aren't currently using KG regions
        #gdf = self.add_kg_attributes(gdf)

        return gdf

    def build_points(self, event_perimeters):
        """Build point geometries of fire events.

        Parameters
        ----------
        event_perimeters : list[firedpy.model_classes.EventPerimeter]
            List of EventPerimeter objects.

        Returns
        -------
        geopandas.geodataframe.GeoDataFrame : A GeoDataFrame of points
            representing grouped fire detections with event IDs, x and y
            coordinates, and detection dates.
        """
        # Build a datetime coordinate dataframe from the events
        logger.info("Converting event perimeters frame to data frame...")
        data = []
        for event in event_perimeters:
            for coords in event.spacetime_coordinates:
                coord0, coord1 = coords[:2]
                coord3 = self._convert_unix_day_to_calendar_date(coords[2])
                data.append([event.event_id, coord1, coord0, coord3])
        df = pd.DataFrame(data, columns=["id", "x", "y", "date"])

        # Center pixel coordinates
        df.loc[:, "x"] = df["x"] + (self.geom[1] / 2)
        df.loc[:, "y"] = df["y"] + (self.geom[-1] / 2)

        # Each entry gets a point object from the x and y coordinates.
        logger.info("Converting event data frame to GeoDataFrame...")
        df.loc[:, "geometry"] = df[["x", "y"]].apply(
            lambda x: Point(tuple(x)),
            axis=1
        )
        gdf = gpd.GeoDataFrame(df, crs=self.crs.proj4, geometry="geometry")

        # Clip to study area if requested
        if self.country:
            shape_file = get_country_file(self.country)
        if shape_file is not None:
            shp_name = Path(shape_file).name
            logger.info(f"Clipping event GeoDataFrame with {shp_name}...")
            gdf = self.clip_to_shape_file(gdf, shape_file)

        return gdf

    def classify_events(self):
        """Classify wildfire events perimeters by tile and merge together."""
        logger.info("Building fire event perimeters.")

        # Build & Collect fire event perimeters
        fire_events = []
        if self.n_cores == 1:
            # Run serially
            for nc_fpath in tqdm(self.files):
                out = self.process_perimeters(nc_fpath)
                fire_events.extend(out)
        else:
            # Run in parallel
            with ProcessPoolExecutor(self.n_cores) as pool:
                jobs = []
                for nc_fpath in self.files:
                    jobs.append(pool.submit(self.process_perimeters, nc_fpath))
                for job in tqdm(as_completed(jobs), total=len(jobs)):
                    fire_events.extend(job.result())

        # Assign event IDs
        for i in range(len(fire_events)):
            fire_events[i].event_id = i

        # Find and merge the edge cases
        if fire_events:
            non_edge = [e for e in fire_events if not e.is_edge]
            edge_events = [e for e in fire_events if e.is_edge]
            for edge_event in edge_events:
                edge_event.compute_min_max()
            merged_edges = self.merge_fire_edge_events(edge_events)
            del edge_events
            last_non_edge_id = non_edge[-1].event_id + 1
            for edge in merged_edges:
                edge.event_id = last_non_edge_id
                last_non_edge_id += 1

            # Combine edge and non edge cases
            fire_events = non_edge + merged_edges

        return fire_events

    def clip_to_shape_file(self, gdf, shape_file):
        """Clip fire events to shapefile boundry, keep overlapping events.

        Parameters
        ----------
        gdf : geopandas.geodataframe.GeoDataFrame
            GeoDataFrame of Firedpy fire events.
        shape_file : str | pathlib.PosixPath
            Path to shapefile representing a study area.

        Returns
        -------
        geopandas.geodataframe.GeoDataFrame : A clipped GeoDataFrame.
        """
        # Read in the shapefile
        shp = gpd.read_file(shape_file)
        shp.to_crs(gdf.crs, inplace=True)

        # Characterize shapefile intersecting events by ID
        shp.loc[:, "intersects"] = 1
        gdf = gpd.sjoin(gdf, shp, how="left")
        gdf.loc[:, "keep"] = gdf.groupby("id")["intersects"].transform("any")

        # Drop non-intersecting events
        clipped_gdf = gdf[gdf["keep"]]
        for tmp_field in ["index_right", "intersects", "keep"]:
            del clipped_gdf[tmp_field]

        return clipped_gdf

    @staticmethod
    def _create_did_column(df, columns):
        """Hashes multiple columns in a DataFrame"""
        df["temp"] = df[columns].apply(
            lambda row: '_'.join(row.values.astype(str)),
            axis=1
        )
        df["did"] = df["temp"].apply(
            lambda x: hashlib.md5(x.encode()).hexdigest()
        )
        df.drop('temp', axis=1, inplace=True)
        return df

    def _create_event_grid_array(self, events: List[EventPerimeter]):
        # Assuming `instances` is your list of class instances and each
        # instance has a `.spacetime_coordinates` attribute that is a list of
        # (y, x, t) tuples.

        # Step 1: Determine min and max coordinates to define the array size.
        min_x = min_y = float('inf')
        max_x = max_y = float('-inf')

        for instance in events:
            for y, x, t in instance.spacetime_coordinates:
                min_x = min(min_x, x)
                max_x = max(max_x, x)
                min_y = min(min_y, y)
                max_y = max(max_y, y)

        min_x_event = [e for e in events if e.min_x == min_x][0]
        max_y_event = [e for e in events if e.max_y == max_y][0]

        max_y_geom = max_y_event.max_geom_y + self._res

        width = min_x_event.min_geom_x
        width_index = 0
        while width < max_x:
            width += self._res
            width_index += 1

        height = max_y_geom
        height_index = 0
        while height > min_y:
            height -= self._res
            height_index += 1

        array_shape = (len(events), height_index+1, width_index+1)
        three_d_array = np.zeros(array_shape, dtype=np.uint16)

        minx = min_x_event.min_geom_x
        miny = max_y_event.max_geom_y
        xs = [math.ceil(minx + self._res * i) for i in range(width_index + 1)]
        ys = [math.ceil(miny - self._res * i) for i in range(height_index + 1)]
        coordinates = {
            'x': np.array(xs),
            'y': np.array(ys)
        }

        for instance_num, instance in enumerate(events):
            for y, x, t in instance.spacetime_coordinates:
                x_idx = int((x - min_x_event.min_geom_x + 1) / self._res)
                y_idx = int((max_y_geom - y - 1) / self._res)
                three_d_array[instance_num, y_idx, x_idx] = t

        return three_d_array, coordinates

    def _extract_date_parts(self, path):
        fname = path.name
        match = re.match(self._file_regex, fname)
        if match:
            match_year = int(match.groupdict()["year"])
            match_day = int(match.groupdict()["ordinal_day"])
            return match_year, match_day
        return None

    @property
    def files(self):
        """Return list of files for given tiles and years."""
        files = [self._generate_local_nc_path(t) for t in self.tiles]
        files = sorted(files)
        if not files:
            raise FileNotFoundError(
                f"Could not find any netcdf files in {self._nc_dir} for "
                f"tiles: {self.tiles}"
            )
        return files

    def get_output_paths(self, project_name, project_directory,
                         start_year, end_year, shape_type):
        """Get dictionary of all output paths for a firedpy run.

        Parameters
        ----------
        project_name : str | NoneType
            A name used to identify the output files of this project.
        project_directory : str
            Project output directory path.
        start_year : int
            The first year of fire events.
        end_year : int
            The last year of fire events.
        shape_type : str
            Build shapefiles from the event data frame. Specify either "shp",
            "gpkg", or both. Shapefiles of both daily progression and overall
            event perimeters will be written to the 'outputs/shapefiles'
            folder of the chosen project directory. These will be saved in the
            specified geopackage format (.gpkg), ERSI Shapefile format (.shp),
            or save them in both formats using the file basename of the fire
            event data frame (e.g. 'modis_events_daily.gpkg' and
            'modis_events.gpkg').

        Returns
        -------
        dict : Dictionary of paths for final firedpy outputs.
        """
        # Set base name
        base_file_name = f"fired_{project_name}_{start_year}_to_{end_year}"

        # Set output directories
        model_outputs_dir = Path(project_directory).joinpath("outputs")
        shape_dir = model_outputs_dir.joinpath("shapefiles")
        table_dir = model_outputs_dir.joinpath("tables")

        # Make output directories
        shape_dir.mkdir(parents=True, exist_ok=True)
        table_dir.mkdir(exist_ok=True)

        # Event-level output paths
        event_base = f"{base_file_name}_events"
        event_csv_path = table_dir.joinpath(f"{event_base}.csv")
        event_shape_path, event_gpkg_path = generate_path(
            project_directory=project_directory,
            base_filename=f"{base_file_name}_events",
            shape_type=shape_type
        )

        # Daily-level output paths
        daily_base = f"{base_file_name}_daily"
        daily_csv_path = table_dir.joinpath(f"{daily_base}.csv")
        daily_shape_path, daily_gpkg_path = generate_path(
            project_directory=project_directory,
            base_filename=daily_base,
            shape_type=shape_type
        )

        # Package output paths
        out_paths = dict(
            event_csv_path=event_csv_path,
            event_shape_path=event_shape_path,
            event_gpkg_path=event_gpkg_path,
            daily_csv_path=daily_csv_path,
            daily_shape_path=daily_shape_path,
            daily_gpkg_path=daily_gpkg_path
        )

        return out_paths

    def group_by_t(self, events) -> List[List[EventPerimeter]]:
        t_groups = []
        for event_group in events:
            events_sorted_by_t = sorted(event_group, key=lambda e: e.min_t)
            t_group = [[events_sorted_by_t[0]]]
            last_group_max_t = t_group[-1][-1].max_t

            for event in events_sorted_by_t[1:]:
                # Handle overlaps in t with the last group
                if event.min_t <= last_group_max_t + self.temporal_param:
                    t_group[-1].append(event)
                    last_group_max_t = max(last_group_max_t, event.max_t)
                else:
                    t_group.append([event])
                    last_group_max_t = event.max_t
            t_groups.extend(t_group)

        return t_groups

    def group_by_x(self, events) -> List[List[EventPerimeter]]:
        events_sorted_by_x = sorted(events, key=lambda event: event.min_x)

        x_groups = [[events_sorted_by_x[0]]]
        last_group_max_x = x_groups[-1][-1].max_x

        for event in events_sorted_by_x[1:]:
            # Overlaps in x with the last group
            if event.min_x <= last_group_max_x + self.sp_buf:
                x_groups[-1].append(event)
                last_group_max_x = max(last_group_max_x, event.max_x)
            else:
                x_groups.append([event])
                last_group_max_x = event.max_x

        return x_groups

    def group_by_y(self, x_groups) -> List[List[EventPerimeter]]:
        y_groups = []
        for group in x_groups:
            group_sorted_by_y = sorted(group, key=lambda event: event.max_y)
            y_group = [[group_sorted_by_y[0]]]
            last_group_min_y = y_group[-1][-1].min_y

            for event in group_sorted_by_y[1:]:
                # Overlaps in y with the last group
                if event.max_y >= last_group_min_y - self.sp_buf:
                    y_group[-1].append(event)
                    last_group_min_y = min(last_group_min_y, event.min_y)
                else:
                    y_group.append([event])
                    last_group_min_y = event.min_y

            y_groups.extend(y_group)

        return y_groups

    @staticmethod
    def _max_growth_date(x: gpd.GeoDataFrame):
        dates = x["date"].to_numpy()
        pixels = x["pixels"].to_numpy()
        loc = np.where(pixels == np.max(pixels))[0]
        d = np.unique(dates[loc])[0]
        return d

    def merge_fire_edge_events(
            self,
            edge_events: List[EventPerimeter]
    ) -> List[EventPerimeter]:
        # Group by close in x and y
        if not edge_events:
            return []

        groups = self.group_by_t(self.group_by_y(self.group_by_x(edge_events)))

        merged_events = []
        logger.info(f"Merging edge {len(groups)} tiles")
        for group in tqdm(groups):
            group_array, coordinates = self._create_event_grid_array(group)
            event_grid = EventGrid(
                project_directory=self.project_directory,
                input_array=group_array,
                coordinates=coordinates
            )
            perimeters = event_grid.get_event_perimeters(all_t=True)
            del event_grid, group_array, coordinates, group
            merged_events.extend(perimeters)

        return merged_events

    def _modis_to_lat_lon(self, x_sinu, y_sinu):
        # Define the MODIS sinusoidal projection (SR-ORG:6974)
        modis_sinu_proj = Proj(
            "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext"
        )

        # Define the WGS84 projection (EPSG:4326)
        wgs84_proj = Proj(proj="latlong", datum="WGS84")

        # Convert from MODIS sinusoidal to WGS84
        lon, lat = transform(modis_sinu_proj, wgs84_proj, x_sinu, y_sinu)

        logger.info(f"Longitude: {lon}, Latitude: {lat}")

    def process_daily_data(self, gdf):
        """Process daily data geometries.

        Parameters
        ----------
        gdf : geopandas.geodataframe.GeoDataFrame
            A firedpy geodataframe of burn events.

        Returns
        -------
        geopandas.geodataframe.GeoDataFrame : The processed geodataframe.
        """
        logger.info("Dissolving polygons...")
        gdfc = gdf.copy()
        gdfc = self._create_did_column(gdfc, ["date", "id"])
        gdfd = gdfc.dissolve(by="did", as_index=False)

        logger.info("Converting polygons to multipolygons...")
        gdfd["geometry"] = gdfd["geometry"].apply(self._as_multi_polygon)

        return gdf

    def process_event_data(self, gdf):
        """Process a firedpy geodataframe for non-daily ouputs.

        Parameters
        ----------
        gdf : geopandas.geodataframe.GeoDataFrame
            A firedpy geodataframe of burn events.

        Returns
        -------
        geopandas.geodataframe.GeoDataFrame : The processed geodataframe.
        """
        edf = gdf.copy()
        edf = edf.drop(["pixels", "date", "event_day", "dy_ar_km2"], axis=1)
        if "did" in gdf.columns:
            edf = edf.drop(["did"])
        if "fid" in edf.columns:
            edf = edf.drop(["fid"])

        logger.info("Dissolving polygons...")
        edf = edf.dissolve(by="id", as_index=False)
        edf.loc[:, "date"] = edf["ig_date"]  # We need the date field later

        logger.info("Calculating perimeter lengths...")
        edf["tot_perim"] = edf["geometry"].to_crs(MODIS_CRS).length

        logger.info("Converting polygons to multipolygons...")
        edf["geometry"] = edf["geometry"].apply(self._as_multi_polygon)

        return edf

    def process_perimeters(self, nc_fpath):
        """Process a perimeter for a given burn data file.

        Parameters
        ----------
        nc_fpath : str
            Path to NetCDF file containing burn data for a MODIS tile.

        Returns
        -------
        list[firedpy.model_classes.EventPerimeter] : A list of fire event
            perimeter objects.
        """
        event_grid = EventGrid(
            nc_fpath=nc_fpath,
            project_directory=self.project_directory,
            spatial_param=self.spatial_param,
            temporal_param=self.temporal_param,
            start_year=self.start_year,
            end_year=self.end_year,
            country=self.country,
            shape_file=self.shape_file
        )
        fire_events = event_grid.get_event_perimeters()
        return fire_events

    def process_geometry(self, gdf):
        """Process event geometries.

        Parameters
        ----------
        gdf : geopandas.geodataframe.GeoDataFrame
            A firedpy geodataframe of burn events.

        Returns
        -------
        geopandas.geodataframe.GeoDataFrame : The geodataframe with processed
            geometries.
        """
        buffer_m = 1 + (self._res / 2)
        logger.info(f"Converting event point to pixel with {buffer_m:.2f}m "
                    "buffer...")
        geometry = gdf.buffer(buffer_m)
        gdf["geometry"] = geometry
        gdf["geometry"] = gdf.envelope
        return gdf

    def save_data(
        self,
        gdf,
        project_name,
        project_directory,
        start_year,
        end_year,
        daily,
        shape_type,
        full_csv
    ):
        """Save event data to various output formats.

        Parameters
        ----------
        gdf : geopandas.geodataframe.GeoDataFrame
            A fully processed firedpy event data frame.
        project_name : str | NoneType
            A name used to identify the output files of this project.
        project_directory : str
            Project output directory path.
        start_year : int
            The first year of fire events.
        end_year : int
            The last year of fire events.
        daily : boolean
            Create the daily polygons or just the event-level perimeter for
            the analysis area. If this flag is set, the daily and event
            polygons will be created, otherwise only the event level.
        shape_type : str
            Build shapefiles from the event data frame. Specify either "shp",
            "gpkg", or both. Shapefiles of both daily progression and overall
            event perimeters will be written to the 'outputs/shapefiles'
            folder of the chosen project directory. These will be saved in the
            specified geopackage format (.gpkg), ERSI Shapefile format (.shp),
            or save them in both formats using the file basename of the fire
            event data frame (e.g. 'modis_events_daily.gpkg' and
            'modis_events.gpkg').
        full_csv : bool
            Write all attributes from input geodataframe to a CSV. Defaults to
            False and includes only a subset of attributes: "x", "y", "id",
            "ig_date", and "last_date".
        """
        # Get all the output file paths
        paths = self.get_output_paths(
            project_name=project_name,
            project_directory=project_directory,
            start_year=start_year,
            end_year=end_year,
            shape_type=shape_type
        )

        # Process and write daily-level events to file if requested
        if daily:
            # Apply output processing for daily version
            ddf = self.process_daily_data(gdf)

            # GeoDataFrames
            if paths["daily_gpkg_path"]:
                dst = paths["daily_gpkg_path"]
                logger.info(f"Saving daily geodataframe to {dst}")
                ddf.to_file(dst, driver="GPKG")
            if paths["daily_shape_path"]:
                dst = paths["daily_shape_path"]
                logger.info(f"Saving daily ESRI Shapefile to {dst}")
                sdf = self.adjust_for_esri(ddf)
                sdf.to_file(dst)

            # CSVs
            dst = paths["daily_csv_path"]
            if full_csv:
                logger.info(f"Writing full daily CSV file to {dst}")
                df = ddf.copy()
                del df["geometry"]
                ddf.to_csv(dst, index=False)
            else:
                logger.info(f"Writing daily CSV file to {dst}")
                df = ddf[["x", "y", "id", "ig_date", "last_date"]]
                df.to_csv(dst, index=False)

        # Process and write event-level events to file
        edf = self.process_event_data(gdf)

        # GeoDataFrames
        if paths["event_gpkg_path"]:
            dst = paths["event_gpkg_path"]
            logger.info(f"Saving event-level geo to {dst}")
            edf.to_file(dst, driver="GPKG")
        if paths["event_shape_path"]:
            dst = paths["event_shape_path"]
            logger.info(f"Writing ESRI Shapefile to {dst}")
            sdf = self.adjust_for_esri(edf)
            sdf.to_file(dst)

        # CSV
        dst = paths["event_csv_path"]
        if full_csv:
            logger.info(f"Writing full daily CSV file to {dst}")
            del edf["geometry"]
            edf.to_csv(dst, index=False)
        else:
            logger.info(f"Writing daily CSV file to {dst}")
            df = edf[["x", "y", "id", "ig_date", "last_date"]]
            df.to_csv(dst, index=False)
