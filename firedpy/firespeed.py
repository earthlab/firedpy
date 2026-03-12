
import math
from logging import getLogger

import numpy as np
import pandas as pd
import shapely
from shapely import LineString
from shapely.geometry import Polygon, MultiPolygon, LineString, Point
from shapely.ops import unary_union, transform, nearest_points
import geopandas as gpd
import pyproj

logger = getLogger(__name__)

def build_cumulative_perims(gdf, id_col="id", date_col="date"):
    gdf = gdf.sort_values([id_col, date_col]).copy()
    gdf[date_col] = pd.to_datetime(gdf[date_col])

    cumulative_list = []

    for fire_id, sub_gdf in gdf.groupby(id_col):
        union_so_far = None
        for geom in sub_gdf.geometry:
            # fix invalid geometries
            geom = geom.buffer(0) if not geom.is_valid else geom

            # cumulative union
            union_so_far = geom if union_so_far is None else shapely.ops.unary_union([union_so_far, geom])
            
            # --- strip all holes robustly ---
            def remove_holes(g):
                if isinstance(g, Polygon):
                    return Polygon(g.exterior)
                elif isinstance(g, MultiPolygon):
                    return MultiPolygon([Polygon(p.exterior) for p in g.geoms])
                else:
                    return g

            union_so_far = remove_holes(union_so_far)

            cumulative_list.append(union_so_far)

    gdf["cum_geom"] = cumulative_list
    return gdf


def computefirespeed(fire_gdf, id_col="id"):
    if fire_gdf.crs is None or fire_gdf.crs.is_geographic:
        raise ValueError(
            f"computefirespeed requires projected CRS in meters, got {fire_gdf.crs}"
        )
    transformer = pyproj.Transformer.from_crs(fire_gdf.crs, "EPSG:4326", always_xy=True)
    geod = pyproj.Geod(ellps="WGS84")
    fire_gdf = fire_gdf.reset_index(drop=True).copy()
    n_rows = fire_gdf.shape[0]
    orig_x = [np.nan] * n_rows
    orig_y = [np.nan] * n_rows
    dest_x = [np.nan] * n_rows
    dest_y = [np.nan] * n_rows
    result_max_dist = [np.nan] * n_rows
    result_speed = [np.nan] * n_rows

    has_ids = id_col in fire_gdf.columns

    ### iterate over time steps
    for i in range(1, fire_gdf.shape[0]):
        # first perimeter in each fire has no valid predecessor
        if has_ids and fire_gdf.iloc[i][id_col] != fire_gdf.iloc[i - 1][id_col]:
            continue

        prev_geom = fire_gdf.iloc[i - 1].cum_geom
        curr_geom = fire_gdf.iloc[i].cum_geom

        #print("timestep:", i)

        # ensure MultiPolygon
        if isinstance(prev_geom, Polygon):
            prev_geom = MultiPolygon([prev_geom])
        if isinstance(curr_geom, Polygon):
            curr_geom = MultiPolygon([curr_geom])

        # --- parent-child intersection matrix ---
        inter_matrix = np.zeros((len(prev_geom.geoms), len(curr_geom.geoms)),dtype=bool)

        for ii in range(len(prev_geom.geoms)):
            for jj in range(len(curr_geom.geoms)):
                inter_matrix[ii, jj] = (
                    prev_geom.geoms[ii].buffer(1e-6).intersects(curr_geom.geoms[jj]))
            
        # --- parent perimeter coordinates ---
        prev_coords = [
            prev_geom.geoms[ii].simplify(0.05).exterior.coords
            for ii in range(len(prev_geom.geoms))
        ]

        best_dist = -np.inf
        best_origin = None
        best_dest = None

        diag_rows = []   # diagnostics for this timestep
        best_child = None

        # --- LOOP OVER CHILD POLYGONS ---
        for j, child_poly in enumerate(curr_geom.geoms):

            parent_ids = np.where(inter_matrix[:, j])[0].tolist()            
            chosen_parent = None

            # --- spot fire handling ---
            if len(parent_ids) == 0:
                dists = [
                    prev_poly.distance(child_poly)
                    for prev_poly in prev_geom.geoms
                ]
                parent_ids = [int(np.argmin(dists))]

            parent_geoms = [prev_geom.geoms[ii] for ii in parent_ids]
            parent_coords = [prev_coords[ii] for ii in parent_ids]

            dist, origin, dest, parent_local_idx = compute_max_vector(
                perim_inner_geoms=parent_geoms,
                perim_outer_geoms=[child_poly],
                inter_matrix=np.ones((len(parent_geoms), 1)), # fix shape: N x 1
                spot_threshold=4000
            )
            
            # infer which parent geometry produced the origin
            chosen_parent = parent_ids[parent_local_idx]

            if dist > best_dist:
                best_dist = dist
                best_origin = origin
                best_dest = dest
                best_child = j

        # --- finalize timestep ---
        if best_origin is None:
            continue
        
        orig_x[i] = best_origin[0]
        orig_y[i] = best_origin[1]
        dest_x[i] = best_dest[0]
        dest_y[i] = best_dest[1]

        lons, lats = transformer.transform(
            [best_origin[0], best_dest[0]],
            [best_origin[1], best_dest[1]]
        )

        dist_m = geod.line_length(lons, lats)

        result_max_dist[i] = dist_m / 1000
        result_speed[i] = (dist_m / 1000) / 24

    return (orig_x, orig_y, dest_x, dest_y, result_max_dist, result_speed)


def compute_max_vector(perim_inner_geoms,
                               perim_outer_geoms,
                               inter_matrix,
                               spot_threshold=4000,
                               debug=False):
    
    result_dist = []
    result_coord_pair = []
    result_poly_pair = []
    result_parent_idx = []
    
    points_per_meter = 1 / 200

    for poly_outer_idx, outer_poly in enumerate(perim_outer_geoms):
        outer_poly = outer_poly.buffer(0)
        spot_flag = not np.any(inter_matrix[:, poly_outer_idx])

        # --------------------------------------------------
        # Determine valid child polygons
        # --------------------------------------------------
        if spot_flag:
            # No intersecting children → compute distances to all children
            distances = [g.distance(outer_poly) for g in perim_inner_geoms]
            nearest_idx = np.argmin(distances)
            if distances[nearest_idx] > spot_threshold:
                continue
            polyids = [nearest_idx]
        else:
            polyids = [ii for ii in range(len(perim_inner_geoms))
                       if inter_matrix[ii, poly_outer_idx]]

        if not polyids:
            continue

        poly_best_dist = -np.inf
        poly_best_pair = None
        poly_best_poly = None
        poly_best_parent_idx = None

        # --------------------------------------------------
        # Iterate over child polygons
        # --------------------------------------------------
        for poly_inner_idx in polyids:
            child_poly = perim_inner_geoms[poly_inner_idx].buffer(0)

            if child_poly.is_empty:
                continue

            # Sample points along child perimeter
            n_child = max(1, int(child_poly.length * points_per_meter))
            if n_child == 0:
                # fallback: use coords from the polygon
                child_pts = [Point(c) for c in child_poly.exterior.coords]
            else:
                # sample_perimeter should return shapely Points
                child_pts = sample_perimeter(child_poly, n_child)

            if child_poly.intersects(outer_poly):
                # overlapping child → compute max distance from child points to parent exterior
                outer_boundary = outer_poly.exterior
                dists = [pt.distance(outer_boundary) for pt in child_pts]
                best_idx = np.argmax(dists)
                pt_child = np.array(child_pts[best_idx].coords[0])
                pt_parent = np.array(outer_boundary.interpolate(outer_boundary.project(child_pts[best_idx])).coords[0])
                max_dist = dists[best_idx]

            else:
                # disconnected child → nearest points as usual
                pt_child_sh, pt_parent_sh = nearest_points(child_poly, outer_poly)
                pt_child = np.array([pt_child_sh.x, pt_child_sh.y])
                pt_parent = np.array([pt_parent_sh.x, pt_parent_sh.y])
                max_dist = np.linalg.norm(pt_parent - pt_child)

            if debug:
                logger.info(f"Poly_outer {poly_outer_idx}, Poly_inner {poly_inner_idx}, "
                    f"dist_val={max_dist:.2f}")
                logger.info(f"Child coord: {pt_child}, Parent coord: {pt_parent}")

            if max_dist > poly_best_dist:
                poly_best_dist = max_dist
                poly_best_pair = (pt_child, pt_parent)
                poly_best_poly = (child_poly, outer_poly)
                poly_best_parent_idx = poly_inner_idx

        # Append results if a valid vector was found
        if poly_best_pair is not None:
            result_dist.append(poly_best_dist)
            result_coord_pair.append(poly_best_pair)
            result_poly_pair.append(poly_best_poly)
            result_parent_idx.append(poly_best_parent_idx)

    # Global maximum across all parents
    if result_dist:
        max_loc = np.argmax(result_dist)
        return (
            result_dist[max_loc],
            result_coord_pair[max_loc][0],
            result_coord_pair[max_loc][1],
            result_parent_idx[max_loc]
        )
    else:
        return np.nan, None, None, None

def sample_perimeter(poly, n_points):
    length = poly.length
    if n_points <= 0:
        return []

    distances = np.linspace(0, length, n_points, endpoint=False)
    sampled_pts = [poly.exterior.interpolate(d) for d in distances]

    # Return as Shapely Points, not NumPy arrays
    return [Point(p.x, p.y) for p in sampled_pts]