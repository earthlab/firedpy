import gc
import os
import re
import shutil
import sys
from collections import OrderedDict, deque
from glob import glob
from typing import List, Set, Dict, Tuple

import numpy as np
import xarray as xr
from shapely import Point, Polygon, MultiPolygon
from tqdm import tqdm
import geopandas as gpd
from datetime import datetime, timedelta
from scipy.spatial import cKDTree
import pandas as pd
import rasterio

from firedpy.__main__ import LandCoverType, EcoRegionType
from firedpy.data_classes import Base

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))


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


class EventPerimeter:
    def __init__(self, event_id, spacetime_coordinates: Set[SpacetimeCoordinate], is_edge: bool = False):
        self.event_id = event_id
        self.merge_id = np.nan
        self.spacetime_coordinates = spacetime_coordinates
        self.is_edge = is_edge

    def add_spacetime_coordinates(self, new_coordinates: Set[SpacetimeCoordinate]):
        self.spacetime_coordinates.update(new_coordinates)

    def __add__(self, other):
        if not isinstance(other, EventPerimeter):
            raise TypeError('Can only add other EventPerimeter')
        self.spacetime_coordinates.update(other.spacetime_coordinates)
        return self


class EventGrid(Base):
    """
    For a single file, however large, find sites with any burn detections
    in the study period, loop through these sites and group by the space-time
    window, save grouped events to a data frame on disk.

    Args:
        spatial_param (int): The allowed distance from the center of the window in index space

    """

    def __init__(self, out_dir: str, nc_file_path: str, spatial_param: int = 5, temporal_param: int = 11,
                 area_unit: str = "Unknown", time_unit: str = "days since 1970-01-01"):
        super().__init__(out_dir)
        self.spatial_param = spatial_param
        self.temporal_param = temporal_param
        self._area_unit = area_unit
        self._time_unit = time_unit
        self._event_grid = {}
        self._next_event_id = 1

        burns = xr.open_dataset(nc_file_path)
        self.data_set = burns
        self._coordinates = burns.coords
        self._input_array = burns.value

    def _add_event_grid(self, event_id: int, new_pts: set):
        self._event_grid.update({p: event_id for p in new_pts})

    def _merge_perimeters(self, perimeters, event_id, obsolete_id):
        # set the merge id in the obsolete id
        perimeters[obsolete_id - 1].merge_id = event_id
        new_pts = perimeters[obsolete_id - 1].spacetime_coordinates

        # update the event_grid and add points to event_id perimeter
        for p in new_pts:
            self._event_grid[p] = event_id
        perimeters[event_id - 1].add_coordinates(new_pts)
        perimeters[event_id - 1].is_edge = perimeters[event_id - 1].is_edge or perimeters[obsolete_id - 1].is_edge

        # set old perimeter to null
        merge_notice = "Merged with event {}".format(event_id)
        perimeters[obsolete_id - 1].spacetime_coordinates = [merge_notice, new_pts]

        return perimeters

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
                origin = coord - param

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
        mask = self.data_set.max(dim="time")
        locs = np.where(mask > 0)
        del mask

        return list(zip(locs[0], locs[1]))

    def get_event_perimeters(self):
        """
        Iterate through each cell in the 3D MODIS Burn Date tile and group it
        into fire events using the space-time window.
        """
        print("Filtering out cells with no events...")
        available_pairs = self._get_available_cells()

        # This is to check the window positions
        nz, ny, nx = self.data_set.shape
        dims = [ny, nx]
        perimeters = []

        # traverse spatially, processing each burn day
        print("Building event perimeters...")
        for pair in tqdm(available_pairs, position=0, file=sys.stdout):
            # Separate coordinates
            y, x = pair

            # get the spatial window
            top, bottom, left, right, center, origin, is_edge = self._get_spatial_window(y, x, dims)
            center_y, center_x = center

            # what if we pull in the window?
            window = self._input_array[:, top:bottom + 1, left:right + 1]

            # The center of the window is the target burn day
            center_burn = window[:, center_y, center_x]
            center_burn = center_burn[center_burn > 0]

            # Loop through each event in the window and identify neighbors
            for burn in center_burn:
                new_pts = set()
                curr_event_ids = list()

                # Now we can get the values and position right away
                val_locs = np.where(abs(burn - window) <= self.temporal_param)
                y_locs = val_locs[1]
                x_locs = val_locs[2]
                oy, ox = origin

                # Now get the geographic coordinates from tile positions
                vals = window[val_locs]
                ys = self._coordinates["y"].data[[oy + yl for yl in y_locs]]
                xs = self._coordinates["x"].data[[ox + xl for xl in x_locs]]

                # Now check if this point is in the event_grid yet
                for i, val in enumerate(vals):
                    curr_pt = SpacetimeCoordinate(ys[i], xs[i], int(val))

                    # already assigned to an event
                    if curr_pt in self._event_grid:
                        curr_event_ids.append(self._event_grid[curr_pt])
                    else:
                        new_pts.add(curr_pt)

                # If this is a new event
                if not curr_event_ids:
                    # create a new perimeter object
                    perimeter = EventPerimeter(self._next_event_id, new_pts, is_edge=is_edge)
                    # append to perimeters list
                    perimeters.append(perimeter)

                    # add points to the grid
                    self._add_event_grid(self._next_event_id, new_pts)

                    # increment the event ID
                    self._next_event_id += 1

                # If all points part of same existing event
                elif len(set(curr_event_ids)) == 1:
                    event_id = curr_event_ids[0]
                    if new_pts:
                        perimeters[event_id - 1].add_spacetime_coordinates(new_pts)
                        self._add_event_grid(event_id, new_pts)

                # events overlap
                else:
                    perimeters = self._merge_perimeters(perimeters,
                                                        curr_event_ids[0],
                                                        curr_event_ids[1])

        return perimeters


class ModelBuilder(Base):
    def __init__(self, out_dir: str, tiles: List[str], spatial_param: int = 5, temporal_param: int = 11):
        super().__init__(out_dir)
        self.tiles = tiles
        self.spatial_param = spatial_param
        self.temporal_param = temporal_param
        self.files = sorted([self._generate_local_nc_path(t) for t in self.tiles])

        if not self.files:
            raise FileNotFoundError(f'Could not find any netcdf files in {self._nc_dir} for tiles: {self.tiles}')

        # Use the first file to get some geometry data for later
        data_set = xr.open_dataset(self.files[0])
        self.crs = data_set.crs
        self.geom = self.crs.geo_transform
        self._res = self.geom[1]
        self.sp_buf = spatial_param * self._res * np.sqrt(2)
        del data_set

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

    def events_within_spatial_range(self, event_1: EventPerimeter, event_2: EventPerimeter) -> bool:
        """
        Check if any point in points_array1 is within distance p to any point in points_array2.

        Args:

        Returns:
            bool: True if any pair of points is closer than p, else False.
        """
        tree1 = cKDTree(np.array([(c.x, c.y) for c in event_1.spacetime_coordinates]))
        tree2 = cKDTree(np.array([(c.x, c.y) for c in event_2.spacetime_coordinates]))

        # query_ball_tree returns indices of all points in tree2 within distance p of points in tree1
        indices = tree1.query_ball_tree(tree2, self.sp_buf, p=2)

        # Check if there are any points within distance p
        for index_list in indices:
            if index_list:  # if not empty, meaning we found a point within the distance
                return True
        return False

    def events_within_temporal_range(self, event_1: EventPerimeter, event_2: EventPerimeter) -> bool:
        event1_burn_days = [c.t for c in event_1.spacetime_coordinates]
        event2_burn_days = np.array([c.t for c in event_2.spacetime_coordinates])

        # EV Find events that are within the start and end of the current edge + temporal param
        return any([(min(event1_burn_days) - self.temporal_param) <= event2_burn_days <=
                    (max(event1_burn_days) + self.temporal_param)])

    def merge_fire_edge_events(self, fire_events: List[EventPerimeter]) -> List[EventPerimeter]:
        # Sort the events by start date
        fire_events = sorted([e for e in fire_events if e.is_edge], key=lambda x: min([c.t for c in
                                                                                       x.spacetime_coordinates]))

        i, j = 0, 1
        events_to_keep = [fire_events[i]]
        while j < len(fire_events):
            if self.events_within_temporal_range(fire_events[i], fire_events[j]):
                if self.events_within_spatial_range(fire_events[i], fire_events[j]):
                    fire_events[i] += fire_events[j]
                else:
                    events_to_keep.append(fire_events[j])
            else:
                events_to_keep.append(fire_events[i])
                i = j
            j += 1

        return events_to_keep

    def build_events(self):
        """
        Use the EventGrid class to classify events tile by tile and then merge
        them all together for a seamless set of wildfire events.

        """
        fire_events: List[EventPerimeter] = []
        for file in self.files:
            fire_events += EventGrid(nc_file_path=file, out_dir=self._out_dir, spatial_param=self.spatial_param,
                                     temporal_param=self.temporal_param).get_event_perimeters()

        fire_events = [e for e in fire_events if not e.is_edge] + self.merge_fire_edge_events(fire_events)

        # Clear memory
        gc.collect()

        return fire_events

    def _clip_to_shape_file(self, gdf: gpd.GeoDataFrame, shape_file_path: str) -> gpd.GeoDataFrame:
        shp = gpd.read_file(shape_file_path)
        shp.to_crs(gdf.crs, inplace=True)
        shp['geometry'] = shp.geometry.buffer(100000)  # Wide buffer to start
        gdf = gpd.sjoin(gdf, shp, how="inner", op="intersects", rsuffix='_join_')
        return gdf.loc[:, ~gdf.columns.str.contains('_join_')]

    def build_points(self, event_perimeter_list: List[EventPerimeter], shape_file_path: str = None):
        # Extract data from EventPerimeter objects and build DataFrame
        df = pd.DataFrame([
            [event.event_id, coord.x, coord.y, self._convert_unix_day_to_calendar_date(coord.t)]
            for event in event_perimeter_list for coord in event.spacetime_coordinates
        ], columns=["id", "x", "y", "date"])

        # Center pixel coordinates
        df["x"] = df["x"] + (self.geom[1] / 2)
        df["y"] = df["y"] + (self.geom[-1] / 2)

        # Each entry gets a point object from the x and y coordinates.
        print("Converting data frame to spatial object...")
        df["geometry"] = df[["x", "y"]].apply(lambda x: Point(tuple(x)), axis=1)
        gdf = gpd.GeoDataFrame(df, crs=self.crs.proj4, geometry=df["geometry"])

        if shape_file_path is not None:
            gdf = self._clip_to_shape_file(gdf, shape_file_path)

        return gdf

    def add_fire_attributes(self, gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        print("Adding fire attributes ...")

        gdf['pixels'] = gdf.groupby(['id', 'date'])['id'].transform('count')

        gdf['ig_utm_x'] = gdf.groupby(['id', 'date'])['x'].transform('nth', 0)

        gdf['ig_utm_y'] = gdf.groupby(['id', 'date'])['y'].transform('nth', 0)

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

        max_date = pd.DataFrame(group[['date', 'pixels']].apply(self._max_growth_date).reset_index())
        max_date = max_date.rename(columns={0: 'mx_grw_dte'})
        gdf = pd.merge(gdf, max_date, on="id")

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

        lc_files = sorted([f for f in glob(os.path.join(self._land_cover_dir, "*tif"))
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
        gdf["did"] = gdf["id"].astype(str) + "-" + gdf["date"]
        return gdf

    def process_daily_data(self, gdf: gpd.GeoDataFrame, output_csv_path: str, daily_shp_path: str,
                           daily_gpkg_path: str) -> gpd.GeoDataFrame:
        print("Dissolving polygons...")
        gdfd = gdf.dissolve(by="did", as_index=False)
        print("Converting polygons to multipolygons...")
        gdfd["geometry"] = gdfd["geometry"].apply(self._as_multi_polygon)

        self.save_data(gdfd, daily_shp_path, daily_gpkg_path, output_csv_path[:-4] + "_daily" + ".csv")

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

        if shape_path is not None:
            gdf.to_file(shape_path)
            print("Saving file to " + shape_path)

        if gpkg_path is not None:
            gdf.to_file(gpkg_path, driver="GPKG")
            print("Saving file to " + gpkg_path)
