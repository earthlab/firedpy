import gc
import math
import os
import re
import shutil
import sys
from collections import OrderedDict, deque
from glob import glob
from itertools import chain
import random
from typing import List, Set, Dict, Tuple
from argparse import Namespace
import multiprocessing as mp
import matplotlib.pyplot as plt
import time

import numpy as np
import xarray as xr
from shapely import Point, Polygon, MultiPolygon
from tqdm import tqdm
import geopandas as gpd
from datetime import datetime
from scipy.spatial import cKDTree
import pandas as pd
import rasterio

from src.enums import LandCoverType, EcoRegionType
from src.data_classes import Base

import cProfile
import pstats

from concurrent.futures import as_completed, ProcessPoolExecutor

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))


def _process_file_perimeter(args):
    event_grid = EventGrid(nc_file_path=args[1], out_dir=args[2],
                           spatial_param=args[4],
                           temporal_param=args[3])
    fire_events = event_grid.get_event_perimeters(args[0], progress_description=os.path.basename(args[1]))
    return fire_events


def _merge_events_spatially_task(namespace: Namespace):
    namespace, pos = namespace
    events: List[EventPerimeter] = namespace.events
    sp_buf: float = namespace.sp_buf

    if len(events) == 1:
        return events

    ckd_trees = [cKDTree(np.array([(c[1], c[0]) for c in event.spacetime_coordinates])) for event in events]

    i, j = 0, 1
    merged_event_indices = set()
    merged_events = []
    pbar = tqdm(total=len(ckd_trees), desc=f"Merging edges for temporal group {pos}", position=pos)
    found_overlapping = False
    while i < len(ckd_trees):
        if j not in merged_event_indices and j < len(ckd_trees):
            indices = ckd_trees[i].query_ball_tree(ckd_trees[j], r=sp_buf, p=2)
            if any(index_list for index_list in indices):
                found_overlapping = True
                events[i].add_spacetime_coordinates(events[j].spacetime_coordinates)
                events[i].is_edge = events[i].is_edge or events[j].is_edge
                ckd_trees[i] = cKDTree(np.array([(c[1], c[0]) for c in events[i].spacetime_coordinates]))
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
        # Ensure t is within 16-bit unsigned integer range, this is ~180 years worth of days so will be OK
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
        # To make the object hashable, return a hash based on the values of x, y, and t
        return hash((self.x, self.y, self.t))

    # The following methods are added for pickling support
    def __getstate__(self):
        return self.x, self.y, self.t

    def __setstate__(self, state):
        self.x, self.y, self.t = state


class EventPerimeter:
    def __init__(self, event_id: int, spacetime_coordinates: Set[Tuple[float, float, np.int16]], is_edge: bool = False):
        self.event_id = event_id
        self.spacetime_coordinates = spacetime_coordinates
        self.is_edge = is_edge
        self.min_y = self.max_y = self.min_x = self.max_x = self.min_t = self.max_t = None
        self.min_geom_x = None
        self.max_geom_y = None

    def add_spacetime_coordinates(self, new_coordinates: Set[Tuple[float, float, np.int16]]):
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
                'min_y': self.min_y, 'max_y': self.max_y,
                'min_x': self.min_x, 'max_x': self.max_x,
                'min_t': self.min_t, 'max_t': self.max_t
            }
        else:
            raise ValueError("Array should be of shape (N, 3) where N is the number of 3-tuples.")

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
        # Hash the tuple of attributes that you're using for equality
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
    """
    For a single file, however large, find sites with any burn detections
    in the study period, loop through these sites and group by the space-time
    window, save grouped events to a data frame on disk.

    Args:
        spatial_param (int): The allowed distance from the center of the window in index space

    """

    def __init__(self, out_dir: str, spatial_param: int = 5, temporal_param: int = 11,
                 area_unit: str = "Unknown", time_unit: str = "days since 1970-01-01", nc_file_path: str = None,
                 input_array: np.array = None, coordinates: Dict[str, np.array] = None):
        super().__init__(out_dir)
        self.spatial_param = spatial_param
        self.temporal_param = temporal_param
        self._area_unit = area_unit
        self._time_unit = time_unit
        self.current_event_id = 0

        if nc_file_path is not None:
            with xr.open_dataset(nc_file_path) as burns:
                self._coordinates = {dim: np.array(burns.coords[dim].values) for dim in ['y', 'x', 'time']}
                self._input_array = burns.value.values
        elif input_array is not None:
            if coordinates is None:
                raise ValueError('If specifying input array, must also provide coordinates parameter')
            self._coordinates = coordinates
            self._input_array = input_array
        else:
            raise ValueError('Must input either nc_file_path or input_array parameter')

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
            """
            Helper function to compute center and origin.
            """
            edges_lower = [t for t in range(param)]
            edges_upper = [dim - t for t in range(param)]

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

        ycenter, oy, is_y_edge = get_edge_positions(array_dims[0], y, self.spatial_param)
        xcenter, ox, is_x_edge = get_edge_positions(array_dims[1], x, self.spatial_param)

        return top, bottom, left, right, [ycenter, xcenter], [oy, ox], is_x_edge or is_y_edge

    def _get_available_cells(self):
        """
        To save time, avoid checking cells with no events at any time step.
        Create a mask of max values at each point. If the maximum at a cell is
        less than or equal to zero there were no values and it will not be
        checked in the event classification step.
        """
        # Using memory - can handle large tiles, but gets pretty high
        locs = np.where(np.any(self._input_array > 0, axis=0))

        return list(zip(locs[0], locs[1]))

    def get_event_perimeters(self, progress_position: int = 0, progress_description: str = '', all_t: bool = False):
        """
        Iterate through each cell in the 3D MODIS Burn Date tile and group it
        into fire events using the space-time window.
        """
        available_pairs = self._get_available_cells()
        nz, ny, nx = self._input_array.shape
        dims = [ny, nx]
        perimeters = RandomAccessSet()
        identified_points = {}

        time_index_buffer = max(1, self.temporal_param // 30)

        for pair in tqdm(available_pairs, position=progress_position, file=sys.stdout, desc=progress_description):
            y, x = pair
            top, bottom, left, right, center, origin, is_edge = self._get_spatial_window(y, x, dims)
            center_y, center_x = center

            window = self._input_array[:, top:bottom + 1, left:right + 1]
            center_burn = window[:, center_y, center_x]

            valid_center_burn_indices = np.where(center_burn > 0)[0]
            valid_center_burn_values = center_burn[valid_center_burn_indices]

            for i, burn_value in enumerate(valid_center_burn_values):
                overlapping_event_ids = RandomAccessSet()
                events_to_remove = set()

                # In each time slice we have a possible range of 30 days (hdf files are monthly temporal resolution)
                l = max(0, int(valid_center_burn_indices[i] - time_index_buffer)) if not all_t else 0
                r = min(nz, int(valid_center_burn_indices[i] + time_index_buffer) + 1) if not all_t else nz

                burn_window = window[l:r, :, :]
                val_locs = np.where(abs(burn_value - burn_window) <= self.temporal_param)
                y_locs, x_locs = val_locs[1], val_locs[2]
                oy, ox = origin

                vals = burn_window[val_locs]

                ys = self._coordinates["y"][oy + y_locs]
                xs = self._coordinates["x"][ox + x_locs]

                new_spacetime_coords = set()
                for i, val in enumerate(vals):
                    c = (np.float32(ys[i]), np.float32(xs[i]), np.uint16(val))

                    if c in identified_points:  # O(1) time
                        overlapping_event_ids.add(identified_points[c])  # O(1) time
                    else:
                        new_spacetime_coords.add(c)

                event = EventPerimeter(event_id=self.current_event_id, spacetime_coordinates=new_spacetime_coords,
                                       is_edge=is_edge)

                if overlapping_event_ids.list:
                    to_keep_id = sorted([e for e in overlapping_event_ids.list], key=lambda l: len(
                        perimeters.get(EventPerimeter(l, set())).spacetime_coordinates))[-1]

                    to_keep = perimeters.get(EventPerimeter(to_keep_id, set()))  # O(1) time
                    for event_id in overlapping_event_ids.set:
                        if event_id != to_keep_id:
                            to_merge = perimeters.get(EventPerimeter(event_id, set()))  # O(1) time
                            to_keep.add_spacetime_coordinates(to_merge.spacetime_coordinates)

                            to_keep.is_edge = to_keep.is_edge or to_merge.is_edge

                            identified_points.update((p, to_keep_id) for p in to_merge.spacetime_coordinates)
                            events_to_remove.add(to_merge)

                    to_keep.add_spacetime_coordinates(event.spacetime_coordinates)
                    to_keep.is_edge = to_keep.is_edge or event.is_edge
                    identified_points.update((p, to_keep_id) for p in new_spacetime_coords)

                    if to_keep.is_edge:
                        if to_keep.min_geom_x is None:
                            to_keep.min_geom_x = self._coordinates['x'][0]
                        if to_keep.max_geom_y is None:
                            to_keep.max_geom_y = self._coordinates['y'][0]

                elif event.spacetime_coordinates:
                    if event.is_edge:
                        event.min_geom_x = self._coordinates['x'][0]
                        event.max_geom_y = self._coordinates['y'][0]

                    perimeters.add(event)  # O(1) time
                    identified_points.update((p, event.event_id) for p in new_spacetime_coords)
                    self.current_event_id += 1

                for event in events_to_remove:
                    perimeters.remove(event)

        return perimeters.list


class ModelBuilder(Base):
    def __init__(self, out_dir: str, tiles: List[str], spatial_param: int = 5, temporal_param: int = 11,
                 n_cores: int = os.cpu_count() - 1):
        super().__init__(out_dir)
        self.tiles = tiles
        self.spatial_param = spatial_param
        self.temporal_param = temporal_param
        self.files = sorted([self._generate_local_nc_path(t) for t in self.tiles])

        if not self.files:
            raise FileNotFoundError(f'Could not find any netcdf files in {self._nc_dir} for tiles: {self.tiles}')

        # Use the first file to get some geometry data for later
        with xr.open_dataset(self.files[0]) as data_set:
            self.crs = data_set.crs
            self.geom = self.crs.geo_transform
            self._res = self.geom[1]
            self.sp_buf = spatial_param * self._res
            self._coordinates = {dim: np.array(data_set.coords[dim].values) for dim in ['y', 'x', 'time']}

        self._n_cores = n_cores

    @staticmethod
    def _max_growth_date(x: gpd.GeoDataFrame):
        dates = x["date"].to_numpy()
        pixels = x["pixels"].to_numpy()
        loc = np.where(pixels == np.max(pixels))[0]
        d = np.unique(dates[loc])[0]
        return d

    @staticmethod
    def _mode(vals) -> float:
        return max(set(list(vals)), key=list(vals).count)

    @staticmethod
    def _as_multi_polygon(polygon):
        if type(polygon) == Polygon:
            polygon = MultiPolygon([polygon])
        return polygon

    @staticmethod
    def _copy_file(src_path: str, dest_path: str) -> str:
        """Generic function to handle copying files."""
        if not os.path.exists(dest_path):
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)
            shutil.copy(src_path, dest_path)
        return dest_path

    def _copy_land_cover_ref(self, land_cover_type: LandCoverType) -> str:
        lookup = os.path.join(PROJECT_DIR, 'ref', 'land_cover', f'MCD12Q1_LegendDesc_Type{land_cover_type.value}.csv')
        land_cover_out_dir_path = os.path.join(self._out_dir, 'tables', 'land_cover',
                                               f"MCD12Q1_LegendDesc_Type{land_cover_type.value}.csv")
        return self._copy_file(lookup, land_cover_out_dir_path)

    def _copy_wwf_file(self) -> str:
        lookup = os.path.join(PROJECT_DIR, 'ref', 'world_eco_regions', 'wwf_terr_ecos.gpkg')
        wwf_out_dir_path = os.path.join(self._out_dir, 'shapefiles', 'eco_region', 'wwf_terr_ecos.gpkg')
        return self._copy_file(lookup, wwf_out_dir_path)

    def _copy_cec_file(self) -> str:
        return self._copy_file(self._project_eco_region_path, self._eco_region_shape_path)

    def _to_kms(self, p: float):
        return (p * self._res ** 2) / 1000000

    def group_by_x(self, events) -> List[List[EventPerimeter]]:
        events_sorted_by_x = sorted(events, key=lambda event: event.min_x)
        x_groups = [[events_sorted_by_x[0]]]

        for event in events_sorted_by_x[1:]:
            if event.min_x <= x_groups[-1][-1].max_x + self.sp_buf:  # Overlaps in x with the last group
                x_groups[-1].append(event)
            else:
                x_groups.append([event])

        return x_groups

    def group_by_t(self, events) -> List[List[EventPerimeter]]:
        t_groups = []
        for group in events:
            group_sorted_by_t = sorted(group, key=lambda event: event.min_t)
            t_group = [[group_sorted_by_t[0]]]

            for event in group_sorted_by_t[1:]:
                if event.min_t <= t_group[-1][-1].max_t + self.temporal_param:  # Overlaps in y with the last group
                    t_group[-1].append(event)
                else:
                    t_group.append([event])

            t_groups.extend(t_group)

        return t_groups

    def group_by_y(self, x_groups) -> List[List[EventPerimeter]]:
        y_groups = []
        for group in x_groups:
            group_sorted_by_y = sorted(group, key=lambda event: event.max_y)
            y_group = [[group_sorted_by_y[0]]]

            for event in group_sorted_by_y[1:]:
                if event.max_y >= y_group[-1][-1].min_y - self.sp_buf:  # Overlaps in y with the last group
                    y_group[-1].append(event)
                else:
                    y_group.append([event])

            y_groups.extend(y_group)

        return y_groups

    def _create_event_grid_array(self, events: List[EventPerimeter]):
        # Assuming `instances` is your list of class instances and each instance
        # has a `.spacetime_coordinates` attribute that is a list of (y, x, t) tuples.

        # Step 1: Determine the min and max coordinates to define the array size.
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

        # Step 2: Convert meter-based coordinates to indices
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

        # Step 3: Create an empty 3D array
        array_shape = (len(events), height_index+1, width_index+1)
        three_d_array = np.zeros(array_shape, dtype=np.uint16)

        coordinates = {
            'x': np.array([math.ceil(min_x_event.min_geom_x + self._res * i) for i in range(width_index + 1)]),
            'y': np.array([math.ceil(max_y_event.max_geom_y - self._res * i) for i in range(height_index + 1)])
        }

        # Step 4: Populate the 3D array with your instances
        for instance_num, instance in enumerate(events):
            for y, x, t in instance.spacetime_coordinates:
                x_idx = int((x - min_x_event.min_geom_x + 1) / self._res)
                y_idx = int((max_y_geom - y - 1) / self._res)
                three_d_array[instance_num, y_idx, x_idx] = t

        return three_d_array, coordinates

    def merge_fire_edge_events(self, edge_events: List[EventPerimeter]) -> List[EventPerimeter]:
        # Group by close in x and y
        if not edge_events:
            return []

        groups = self.group_by_t(self.group_by_y(self.group_by_x(edge_events)))
        #groups = self.group_by_y(self.group_by_x(edge_events))

        merged_events = []
        for i, spatial_group in enumerate(groups):
            print(len(spatial_group), len(set([s.event_id for s in spatial_group])), 'set length')
            group_array, coordinates = self._create_event_grid_array(spatial_group)
            print(sys.getsizeof(group_array), 'bytes', group_array.shape)
            event_grid = EventGrid(out_dir=self._out_dir, input_array=group_array, coordinates=coordinates)
            print(sys.getsizeof(event_grid), 'event grid bytes')
            perimeters = event_grid.get_event_perimeters(progress_position=i,
                                                         progress_description=f'Merging edge tiles for group {i} of'
                                                                              f' {len(groups)}', all_t=True)
            print(sys.getsizeof(perimeters), 'perimeters bytes')
            del event_grid, group_array, coordinates, spatial_group
            merged_events.extend(perimeters)

        print(len(merged_events), 'merged')

        return merged_events

    def build_events(self):
        """
        Use the EventGrid class to classify events tile by tile and then merge
        them all together for a seamless set of wildfire events.

        """
        print('Building fire event perimeters')

        fire_events = []

        # Create a pool of workers with the number of cores specified
        with mp.Pool(self._n_cores) as pool:
            # Map the files to the processing function and collect the results
            # Use starmap to pass the arguments unpacked
            args = []
            for i, file in enumerate(self.files):
                args.append((i, file, self._out_dir, self.temporal_param, self.spatial_param))
            results = pool.imap_unordered(_process_file_perimeter, args)

            # As each file is processed, results will be appended to the fire_events list
            for result in results:
                fire_events.extend(result)

        for i in range(len(fire_events)):
            fire_events[i].event_id = i

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

        fire_events = non_edge + merged_edges

        gc.collect()

        return fire_events

    def _clip_to_shape_file(self, gdf: gpd.GeoDataFrame, shape_file_path: str) -> gpd.GeoDataFrame:
        shp = gpd.read_file(shape_file_path)
        shp.to_crs(gdf.crs, inplace=True)
        shp['geometry'] = shp.geometry.buffer(100000)  # Wide buffer to start
        clipped_gdf = gdf.clip(shp)
        return clipped_gdf

    def build_points(self, event_perimeter_list: List[EventPerimeter], shape_file_path: str = None):
        # Extract data from EventPerimeter objects and build DataFrame
        df = pd.DataFrame([
            [event.event_id, coord[0], coord[1], self._convert_unix_day_to_calendar_date(coord[2])]
            for event in event_perimeter_list for coord in event.spacetime_coordinates
        ], columns=["id", "x", "y", "date"])

        # Center pixel coordinates
        df["x"] = df["x"] + (self.geom[1] / 2)
        df["y"] = df["y"] + (self.geom[-1] / 2)

        # Each entry gets a point object from the x and y coordinates.
        print("Converting data frame to spatial object...")
        df["geometry"] = df[["y", "x"]].apply(lambda x: Point(tuple(x)), axis=1)
        gdf = gpd.GeoDataFrame(df, crs=self.crs.proj4, geometry=df["geometry"])

        if shape_file_path is not None:
            gdf = self._clip_to_shape_file(gdf, shape_file_path)
        return gdf

    def add_fire_attributes(self, gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        print("Adding fire attributes ...")

        # date_id_group = gdf.groupby(['id', 'date'])
        print('finished grouping')

        gdf['pixels'] = gdf.groupby(['id', 'date'])['id'].transform('count')

        gdf['ig_utm_x'] = gdf.groupby(['id', 'date'])['x'].nth(0)

        gdf['ig_utm_y'] = gdf.groupby(['id', 'date'])['y'].nth(0)

        print('grouping')
        group = gdf.groupby('id')

        gdf['date'] = gdf['date'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))

        gdf['ig_date'] = group['date'].transform('min')

        gdf['ig_day'] = gdf['ig_date'].apply(lambda x: datetime.strftime(x, '%j'))

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
        print('to_csv')

        max_date = pd.DataFrame(group[['date', 'pixels']].apply(self._max_growth_date).reset_index())
        max_date = max_date.rename(columns={0: 'mx_grw_dte'})
        gdf = gdf.merge(max_date, on="id")

        print('merged')

        gdf = gdf[['id', 'date', 'ig_date', 'ig_day', 'ig_month',
                   'ig_year', 'last_date', 'event_day', 'event_dur',
                   'pixels', 'tot_pix', 'dy_ar_km2', 'tot_ar_km2',
                   'fsr_px_dy', 'fsr_km2_dy',
                   'mx_grw_px', 'mn_grw_px', 'mu_grw_px',
                   'mx_grw_km2', 'mn_grw_km2', 'mu_grw_km2', 'mx_grw_dte',
                   'x', 'y', 'geometry', 'ig_utm_x', 'ig_utm_y']]

        gdf = gdf.reset_index(drop=True)

        return gdf

    def add_land_cover_attributes(self, gdf: gpd.GeoDataFrame,
                                  land_cover_type: LandCoverType = LandCoverType.NONE) -> gpd.GeoDataFrame:
        print('Adding land cover attributes...')

        # We'll need to specify which type of land_cover
        lc_descriptions = {LandCoverType.IGBP: "IGBP global vegetation classification scheme",
                           LandCoverType.UMD: "University of Maryland (UMD) scheme",
                           LandCoverType.MODIS_LAI: "MODIS-derived LAI/fPAR scheme",
                           LandCoverType.MODIS_NPP: "MODIS-derived Net Primary Production (NPP) scheme",
                           LandCoverType.PFT: "Plant Functional Type (PFT) scheme."}

        lc_files = sorted([f for f in glob(os.path.join(self._land_cover_dir, '**', "*tif"), recursive=True)
                           if re.match(self._land_cover_regex, os.path.basename(f))])

        lc_years = [re.match(self._land_cover_regex, os.path.basename(f)).groupdict()['year'] for f in lc_files]
        lc_files = {lc_years[i]: lc_files[i] for i in range(len(lc_files))}

        # Rasterio point querier (will only work here)
        def point_query(row):
            x = row['x']
            y = row['y']
            try:
                val = int([val for val in lc.sample([(x, y)])][0])
            except Exception:
                val = np.nan
            return val

        # Get the range of burn years
        burn_years = list(gdf['ig_year'].unique())

        # This works faster when split by year and the pointer is outside
        # This is also not the best way
        sgdfs = []
        for year in tqdm(burn_years, position=0, file=sys.stdout):

            sgdf = gdf[gdf['ig_year'] == year]

            # Now set year one back for land_cover
            year = year - 1

            # Use previous year's lc
            if year < min(lc_years):
                year = min(lc_years)
            elif year > max(lc_years):
                year = max(lc_years)

            lc_file = lc_files[year]
            lc = rasterio.open(lc_file)
            sgdf['lc_code'] = sgdf.apply(point_query, axis=1)
            sgdfs.append(sgdf)

        gdf = pd.concat(sgdfs)
        gdf = gdf.reset_index(drop=True)
        gdf['lc_mode'] = gdf.groupby('id')['lc_code'].transform(self._mode)

        # Add in the class description from land_cover tables
        if land_cover_type != LandCoverType.NONE:
            land_cover_path = self._copy_land_cover_ref(land_cover_type)
            lc_table = pd.read_csv(land_cover_path)
            gdf = pd.merge(left=gdf, right=lc_table, how='left', left_on='lc_mode', right_on='Value')
            gdf = gdf.drop('Value', axis=1)
            gdf['lc_type'] = lc_descriptions[land_cover_type]

        gdf.rename({'lc_description': 'lc_desc'}, inplace=True, axis='columns')

        return gdf

    def _add_attributes_from_na_cec(self, gdf, eco_region_level):
        # Different levels have different sources
        eco_types = {
            'NA_L3CODE': ('Level III Ecoregions ' + '(NA-Commission for Environmental Cooperation)'),
            'NA_L2CODE': ('Level II Ecoregions ' + '(NA-Commission for Environmental Cooperation)'),
            'NA_L1CODE': ('Level I Ecoregions ' + '(NA-Commission for Environmental Cooperation)')}

        # Read in the Level File (contains every level) and reference table
        shp_path = self._copy_cec_file()
        eco = gpd.read_file(shp_path)
        eco.to_crs(gdf.crs, inplace=True)

        # Filter for selected level (level III defaults to US-EPA version)
        if not eco_region_level:
            eco_region_level = 3

        eco_code = [c for c in eco.columns if str(eco_region_level) in
                    c and 'CODE' in c]

        if len(eco_code) > 1:
            eco_code = [c for c in eco_code if 'NA' in c][0]
        else:
            eco_code = eco_code[0]

        print("Selected ecoregion code: " + str(eco_code))

        # Find modal eco-region for each event id
        eco = eco[[eco_code, 'geometry']]
        gdf = gpd.sjoin(gdf, eco, how="left", op="within")
        gdf = gdf.reset_index(drop=True)
        gdf['eco_mode'] = gdf.groupby('id')[eco_code].transform(self._mode)

        # Add in the type of eco-region
        gdf['eco_type'] = eco_types[eco_code]

        # Add in the name of the modal ecoregion
        eco_ref = pd.read_csv(self._eco_region_csv_path)
        eco_name = eco_code.replace('CODE', 'NAME')
        eco_df = eco_ref[[eco_code, eco_name]].drop_duplicates()
        eco_map = dict(zip(eco_df[eco_code], eco_df[eco_name]))
        gdf['eco_name'] = gdf[eco_code].map(eco_map)

        # Clean up column names
        gdf = gdf.drop('index_right', axis=1)
        gdf = gdf.drop(eco_code, axis=1)

        return gdf

    def _add_attributes_from_wwf(self, gdf):
        # Read in the world ecoregions from WWF
        eco_path = self._copy_wwf_file()
        eco = gpd.read_file(eco_path)
        eco.to_crs(gdf.crs, inplace=True)

        # Find modal eco region for each event id
        eco = eco[["ECO_NUM", "ECO_NAME", "geometry"]]
        gdf = gpd.sjoin(gdf, eco, how="left", op="within")
        gdf = gdf.reset_index(drop=True)

        gdf["eco_mode"] = gdf.groupby('id')['ECO_NUM'].transform(self._mode)
        gdf["eco_name"] = gdf["ECO_NAME"]
        gdf["eco_type"] = "WWF Terrestrial Ecoregions of the World"

        gdf = gdf.drop('index_right', axis=1)
        gdf = gdf.drop('ECO_NAME', axis=1)
        gdf = gdf.drop('ECO_NUM', axis=1)

        return gdf

    def add_eco_region_attributes(self, gdf: gpd.GeoDataFrame, eco_region_type: EcoRegionType = EcoRegionType.NA,
                                  eco_region_level: int = None) -> gpd.GeoDataFrame:

        print("Adding eco-region attributes ...")
        if eco_region_level or (eco_region_type == EcoRegionType.NA):
            gdf = self._add_attributes_from_na_cec(gdf, eco_region_level)
        else:
            gdf = self._add_attributes_from_wwf(gdf)

        return gdf

    def process_geometry(self, gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        print("Creating buffer...")
        geometry = gdf.buffer(1 + (self._res / 2))
        gdf["geometry"] = geometry
        gdf["geometry"] = gdf.envelope
        gdf["did"] = gdf["id"].astype(str) + "-" + str(gdf["date"])
        return gdf

    def process_daily_data(self, gdf: gpd.GeoDataFrame, output_csv_path: str, daily_shp_path: str,
                           daily_gpkg_path: str) -> gpd.GeoDataFrame:
        print("Dissolving polygons...")
        gdfd = gdf.dissolve(by="did", as_index=False)
        print("Converting polygons to multipolygons...")
        gdfd["geometry"] = gdfd["geometry"].apply(self._as_multi_polygon)

        self.save_data(gdfd, daily_shp_path, daily_gpkg_path, output_csv_path)

        gdf = gdfd.drop(['did', 'pixels', 'date', 'event_day', 'dy_ar_km2'], axis=1)
        gdf = gdf.dissolve(by="id", as_index=False)
        print("Calculating perimeter lengths...")
        gdf["tot_perim"] = gdf["geometry"].length
        gdf["geometry"] = gdf["geometry"].apply(self._as_multi_polygon)

        return gdf

    def process_event_data(self, gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        gdf = gdf.drop(['did', 'pixels', 'date', 'event_day', 'dy_ar_km2'], axis=1)

        print("Dissolving polygons...")
        gdf = gdf.dissolve(by="id", as_index=False)
        print("Calculating perimeter lengths...")
        gdf["tot_perim"] = gdf["geometry"].length
        print("Converting polygons to multipolygons...")
        gdf["geometry"] = gdf["geometry"].apply(self._as_multi_polygon)

        return gdf

    def save_event_data(self, gdf: gpd.GeoDataFrame, output_csv_path: str, event_shape_path: str, event_gpkg_path: str,
                        full_csv: bool):
        event_csv = output_csv_path[:-4] + "_events" + ".csv"
        if full_csv:
            gdf.to_csv(event_csv, index=False)
        else:
            to_raw_csv = gdf[["x", "y", "id", "ig_date", "last_date"]]
            to_raw_csv.to_csv(event_csv, index=False)

        self.save_data(gdf, event_shape_path, event_gpkg_path)

    @staticmethod
    def save_data(gdf: gpd.GeoDataFrame, shape_path: str = None, gpkg_path: str = None, csv_path: str = None):
        if csv_path:
            gdf.to_csv(csv_path, index=False)

        if gpkg_path is not None:
            gdf.to_file(gpkg_path, driver="GPKG")
            print("Saving file to " + gpkg_path)

        if shape_path is not None:
            print('Writing shape file')
            gdf['ig_date'] = [str(d) for d in gdf['ig_date']]
            gdf['last_date'] = [str(d) for d in gdf['last_date']]
            gdf['mx_grw_dte'] = [str(d) for d in gdf['mx_grw_dte']]

            gdf.to_file(shape_path)
            print("Saving file to " + shape_path)


