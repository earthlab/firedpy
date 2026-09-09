
import math
from logging import getLogger

import numpy as np
import pandas as pd
import shapely
from shapely import LineString
from shapely.geometry import Polygon, MultiPolygon, LineString, Point
from shapely.ops import unary_union, transform, nearest_points
from shapely.strtree import STRtree
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


def computefirespeed(fire_gdf,id_col="id",spot_threshold=20_000):
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

    for i in range(1, fire_gdf.shape[0]):
        
        if has_ids and fire_gdf.iloc[i][id_col] != fire_gdf.iloc[i - 1][id_col]:
            continue

        prev_geom = fire_gdf.iloc[i - 1].cum_geom
        curr_geom = fire_gdf.iloc[i].cum_geom

        if isinstance(prev_geom, Polygon):
            prev_geom = MultiPolygon([prev_geom])
        if isinstance(curr_geom, Polygon):
            curr_geom = MultiPolygon([curr_geom])

        inter_matrix = np.zeros((len(prev_geom.geoms), len(curr_geom.geoms)), dtype=bool)
        prev_buffered = [p.buffer(1e-6) for p in prev_geom.geoms]
        for ii in range(len(prev_geom.geoms)):
            for jj in range(len(curr_geom.geoms)):
                inter_matrix[ii, jj] = prev_buffered[ii].intersects(curr_geom.geoms[jj])

        prev_coords = [
            prev_geom.geoms[ii].simplify(0.05).exterior.coords
            for ii in range(len(prev_geom.geoms))
        ]

        # --- Skip if geometry is identical between timesteps ---
        if prev_geom.equals(curr_geom):
            result_max_dist[i] = 0
            result_speed[i] = 0
            # orig_x/y and dest_x/y remain NaN
            continue

        best_dist = -np.inf
        best_origin = None
        best_dest = None
        best_child = None

        for j, child_poly in enumerate(curr_geom.geoms):
            parent_ids = np.where(inter_matrix[:, j])[0].tolist()
            spot = len(parent_ids) == 0
            if spot:
                dists = [prev_poly.distance(child_poly) for prev_poly in prev_geom.geoms]
                parent_ids = [int(np.argmin(dists))]

            parent_geoms = [prev_geom.geoms[ii] for ii in parent_ids]
            parent_coords = [prev_coords[ii] for ii in parent_ids]

            dist, origin, dest, parent_local_idx = compute_max_vector(
                perim_inner_geoms=parent_geoms,
                perim_outer_geoms=[child_poly],
                inter_matrix=inter_matrix[parent_ids][:, [j]],
                # For spots, don't apply the 20-km cutoff in MODIS sinusoidal coordinates.
                spot_threshold=np.inf if spot else spot_threshold,)
            
            # Final SPOT threshold is GEODESIC.
            if spot and origin is not None and dest is not None:
                lons, lats = transformer.transform([origin[0], dest[0]],[origin[1], dest[1]],)
                spot_dist_m = geod.line_length(lons, lats)
                if spot_dist_m > spot_threshold:
                    continue

            if dist > best_dist:
                best_dist = dist
                best_origin = origin
                best_dest = dest
                best_child = j

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
        result_speed[i] = (dist_m / 1000) # in km/day

    return (orig_x, orig_y, dest_x, dest_y, result_max_dist, result_speed)


def compute_max_vector(perim_inner_geoms,
                       perim_outer_geoms,
                       inter_matrix,
                       spot_threshold=20000):

    result_dist = []
    result_coord_pair = []
    result_poly_pair = []
    result_parent_idx = []

    points_per_meter = 1 / 200

    for poly_outer_idx, outer_poly in enumerate(perim_outer_geoms):
        outer_poly = outer_poly.buffer(0)
        spot_flag = not np.any(inter_matrix[:, poly_outer_idx])

        if spot_flag:
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

        # --- Sample child boundary once ---
        n_child_pts = max(1, int(outer_poly.length * points_per_meter))
        child_pts_sample = sample_perimeter(outer_poly, n_child_pts)
        child_coords_arr = np.array([p.coords[0] for p in child_pts_sample])  # (N_child, 2)

        # --- Check if any parent intersects (vs all spots) ---
        any_intersecting = any(
            perim_inner_geoms[ii].intersects(outer_poly) for ii in polyids
        )

        if any_intersecting:
            # -------------------------------------------------------
            # OVERLAPPING CASE: exact nearest-boundary maximin
            # -------------------------------------------------------

            # Build a combined set of densely-sampled parent boundary points
            # and use STRtree for fast nearest-neighbor lookup
            all_parent_pts = []
            all_parent_pt_ids = []
            all_parent_polys_for_check = []
            
            for ii in polyids:
                parent_poly = perim_inner_geoms[ii].buffer(0)
                if parent_poly.is_empty or not parent_poly.intersects(outer_poly):
                    continue
                n_pts = max(1, int(parent_poly.exterior.length * points_per_meter))
                pts = sample_perimeter(parent_poly, n_pts)
                all_parent_pts.extend(pts)
                all_parent_pt_ids.extend([ii] * len(pts))
                all_parent_polys_for_check.extend([parent_poly] * len(pts))
            
            if not all_parent_pts:
                continue
            
            # Build STRtree on parent boundary points
            tree = STRtree(all_parent_pts)

            # Bulk nearest-neighbor query — replaces per-point loop
            nearest_idxs = tree.nearest(child_pts_sample)  # array of indices, one per child pt

            # Pre-compute which child points are inside any parent (need wrong-side check)
            all_parent_union = unary_union([perim_inner_geoms[ii].buffer(0) for ii in polyids])
            child_inside_mask = np.array([all_parent_union.contains(pt) for pt in child_pts_sample])
            child_outside_mask = ~child_inside_mask  # reuse for fallback

            exact_D_min = np.zeros(len(child_pts_sample))
            exact_nearest_parent_pts = []
            exact_nearest_parent_ids = []

            for ci, (pt_c, nearest_idx) in enumerate(zip(child_pts_sample, nearest_idxs)):
                near_pt = all_parent_pts[nearest_idx]
                parent_poly_b = all_parent_polys_for_check[nearest_idx]
                best_id = all_parent_pt_ids[nearest_idx]

                # Wrong-side check only for child points inside a parent
                if child_inside_mask[ci]:
                    test_line = LineString([near_pt.coords[0], pt_c.coords[0]])
                    if test_line.length > 0:
                        interior_overlap = test_line.intersection(parent_poly_b)
                        interior_fraction = (interior_overlap.length / test_line.length
                                            if not interior_overlap.is_empty else 0)
                        if interior_fraction > 0.1:
                            boundary = parent_poly_b.exterior
                            hits = test_line.intersection(boundary)
                            if not hits.is_empty:
                                near_pt = (min(hits.geoms, key=lambda p: p.distance(pt_c))
                                          if hasattr(hits, 'geoms') else hits)

                exact_D_min[ci] = pt_c.distance(near_pt)
                exact_nearest_parent_pts.append(near_pt)
                exact_nearest_parent_ids.append(best_id)

            # maximin: child point whose nearest parent boundary point is farthest
            best_child_idx = np.argmax(exact_D_min)
            pt_child = np.array(child_pts_sample[best_child_idx].coords[0])
            pt_parent = np.array(exact_nearest_parent_pts[best_child_idx].coords[0])
            chosen_parent_id = exact_nearest_parent_ids[best_child_idx]
            max_dist = exact_D_min[best_child_idx]

            # --- Validate: 75% of vector must be inside child ---
            test_line = LineString([tuple(pt_parent), tuple(pt_child)])

            if test_line.length == 0:
                continue

            child_overlap = test_line.intersection(outer_poly)
            child_length = child_overlap.length if not child_overlap.is_empty else 0
            child_fraction = child_length / test_line.length

            if child_fraction >= 0.75:
                poly_best_dist = max_dist
                poly_best_pair = (pt_parent, pt_child)
                poly_best_poly = (perim_inner_geoms[chosen_parent_id].buffer(0), outer_poly)
                poly_best_parent_idx = chosen_parent_id

            else:
                # --- Fallback: try each child point in descending order of D_min,
                #     find first one whose vector passes the 75% overlap check ---

                sorted_child_rows = np.argsort(exact_D_min)[::-1]  # descending by distance
                outside_rows = [ci for ci in sorted_child_rows if child_outside_mask[ci]]
                rows_to_try = outside_rows if outside_rows else sorted_child_rows  # fallback to all if none outside

                found = False
                for ci in rows_to_try:  
                    candidate_pt_child = np.array(child_pts_sample[ci].coords[0])
                    candidate_pt_parent = np.array(exact_nearest_parent_pts[ci].coords[0])
                    candidate_line = LineString([tuple(candidate_pt_parent), tuple(candidate_pt_child)])
                    if candidate_line.length == 0:
                        continue

                    overlap = candidate_line.intersection(outer_poly)
                    overlap_length = overlap.length if not overlap.is_empty else 0
                    if overlap_length / candidate_line.length >= 0.75:
                        pt_parent = candidate_pt_parent
                        pt_child = candidate_pt_child
                        chosen_parent_id = exact_nearest_parent_ids[ci]
                        max_dist = exact_D_min[ci]
                        found = True
                        break

                if not found:
                    # Last resort: use best maximin result regardless of overlap
                    pt_child = np.array(child_pts_sample[best_child_idx].coords[0])
                    pt_parent = np.array(exact_nearest_parent_pts[best_child_idx].coords[0])
                    chosen_parent_id = exact_nearest_parent_ids[best_child_idx]
                    max_dist = exact_D_min[best_child_idx]

                poly_best_dist = max_dist
                poly_best_pair = (pt_parent, pt_child)
                poly_best_poly = (perim_inner_geoms[chosen_parent_id].buffer(0), outer_poly)
                poly_best_parent_idx = chosen_parent_id

        else:
            # -------------------------------------------------------
            # DISCONNECTED / SPOT FIRE CASE
            # -------------------------------------------------------
            for ii in polyids:
                parent_poly = perim_inner_geoms[ii].buffer(0)
                if parent_poly.is_empty:
                    continue

                pt_parent_sh, pt_child_sh = nearest_points(parent_poly, outer_poly)
                parent_anchor = Point(pt_parent_sh.x, pt_parent_sh.y)
                dists = [pt.distance(parent_anchor) for pt in child_pts_sample]
                best_idx = np.argmax(dists)
                pt_parent = np.array([parent_anchor.x, parent_anchor.y])
                pt_child = np.array(child_pts_sample[best_idx].coords[0])
                max_dist = dists[best_idx]

                # Re-apply spot threshold on max-distance (applied to min distance before)
                if max_dist > spot_threshold:
                    continue

                if max_dist > poly_best_dist:
                    poly_best_dist = max_dist
                    poly_best_pair = (pt_parent, pt_child)
                    poly_best_poly = (parent_poly, outer_poly)
                    poly_best_parent_idx = ii

        if poly_best_pair is not None:
            result_dist.append(poly_best_dist)
            result_coord_pair.append(poly_best_pair)
            result_poly_pair.append(poly_best_poly)
            result_parent_idx.append(poly_best_parent_idx)

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
