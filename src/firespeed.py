import math
import numpy as np
from shapely import LineString
import geopandas as gpd
import pyproj


def computefirespeed(fire_gdf):
    transformer = pyproj.Transformer.from_crs(fire_gdf.crs, "EPSG:4326", always_xy=True)
    geod = pyproj.Geod(ellps="WGS84")
    orig_x = [np.nan]
    orig_y = [np.nan]
    dest_x = [np.nan]
    dest_y = [np.nan]
    result_max_dist = [np.nan]
    result_speed = [np.nan]
    prev_step = [fire_gdf.iloc[0]["geometry"].geoms[ii].simplify(0.05).exterior.coords for ii in range(len(fire_gdf.iloc[0]["geometry"].geoms))]
    ### need to deal with resampling or something here in future iteration of code
    ### iterate over time steps
    for i in range(1, min(fire_gdf.shape[0], 10)):
        ### setup for overlap and spot checks...
        inter_matrix = np.zeros((len(fire_gdf.iloc[i-1]["geometry"].geoms), len(fire_gdf.iloc[i]["geometry"].geoms)))
        for ii in range(inter_matrix.shape[0]):
            for jj in range(inter_matrix.shape[1]):
                inter_matrix[ii, jj] = fire_gdf.iloc[i-1]["geometry"].geoms[ii].intersects(fire_gdf.iloc[i]["geometry"].geoms[jj])
        ### bug check
        if inter_matrix.shape[0] == 1 and inter_matrix.shape[1] == 1 and inter_matrix[0, 0] == False:
            print("Error: No overlap at", i)

        ### return maximum_distance, max_dist_origin, max_dist_destination
        max_fire_dist, max_origin, max_destination, prev_step = compute_max_vector(fire_gdf["geometry"][i-1].geoms, fire_gdf["geometry"][i].geoms, prev_step, 
                                                                                   inter_matrix, buffer=200, maxbins=200, slop=2)
        orig_x.append(max_origin[0])
        orig_y.append(max_origin[1])
        dest_x.append(max_destination[0])
        dest_y.append(max_destination[1])

        ### compute meter distance
        lons, lats = transformer.transform([max_origin[0], max_destination[0]],
                                           [max_origin[1], max_destination[1]])
        dist = geod.line_length(lons, lats)
        result_max_dist.append(dist)

        ### assumption here is these are daily perims...
        ### so to compute spread in km/h, we do (dist in km) / 24
        result_speed.append((dist * 1000) / 24)

    
        
    return orig_x, orig_y, dest_x, dest_y, result_max_dist, result_speed


def compute_max_vector(perim_inner_geoms, perim_outer_geoms, inner_coords, inter_matrix, buffer, maxbins, slop):
    ### distance computation 2 -- binned
    outer_coords = []

    ### does computing this beforehand speed things up? maybe a bit
    root2 = round(math.sqrt(2), 5)

    result_dist = []
    result_coord_pair = []
    result_poly_pair = []
    ### perim_inner/outer geoms to be just gpd.iloc[i/i+1].geoms
    ### broadly, for points vi, vj and polys Px Py,
    ### minimize distances between vi and vj with i fixed (closest point to fixed point)
    ### maximize distances between vi and vjmax with vjmax the closest point to vi (furthest travelled between 2 polys)
    ### minimize distances between Px and Py with x fixed (closest polygon to the one looked at)
    ### maximize distances between Px and Pymax (furthest point from closest point on closest polygon)
    for poly_outer in range(len(perim_outer_geoms)):
        buffer_poly = perim_outer_geoms[poly_outer].buffer(buffer)
        ### computer poly_outer bounding box
        outer_bbox = perim_outer_geoms[poly_outer].exterior.bounds
        outer = perim_outer_geoms[poly_outer].simplify(0.05).exterior.coords
        #else:
        #    outer = resample(perim_outer_geoms[poly_outer].simplify(0.05).exterior, params["resample"]).coords
        outer_coords.append(outer)
        ### this checks if any previous perimeter is inside this one...
        ### if not, spotted
        spot_flag = not np.any(inter_matrix[:, poly_outer])
        polyids = []
        if spot_flag:
            ### in the R code, this is resampling to get more points per meter to get accurate estimate...?
            ### need to compare with all polygons
            polyids= [ii for ii in range(len(perim_inner_geoms))]
        ### if it hasn't spotted, compare to all perims that fall inside:
        else:
            for ii in range(len(perim_inner_geoms)):
                if inter_matrix[ii, poly_outer]:
                    polyids.append(ii)
        ### now, keep track of all comparisons for this polygon
        poly_min_dist = float("inf")
        poly_coordpair = None
        poly_pair = None

        for poly_inner in polyids:
            ### compute inner bounding box...
            inner_bbox = perim_inner_geoms[poly_inner].exterior.bounds
            ### compute combined bounding box
            bbox = (min(inner_bbox[0], outer_bbox[0]), min(inner_bbox[1], outer_bbox[1]), 
                    max(inner_bbox[2], outer_bbox[2]), max(inner_bbox[3], outer_bbox[3]))
            ### guess x, y bin sizes 
            bins_x = math.ceil((bbox[2] - bbox[0])/maxbins)
            bins_y = math.ceil((bbox[3] - bbox[1])/maxbins)
            ### bin resolution is bigger of these values since we want square bins
            bin_res = max(bins_x, bins_y)
            ### compute actual number of bins with bin resolution 
            grid_size = (math.ceil((bbox[2] - bbox[0])/bin_res), math.ceil((bbox[3] - bbox[1])/bin_res))
            grid_spatial = (grid_size[0] * bin_res, grid_size[1] * bin_res)
            grid_offset = ((grid_spatial[0] - (bbox[2] - bbox[0]))/2, (grid_spatial[1] - (bbox[3] - bbox[1]))/2)

            ### now make np array w/ coarser side...
            inner_bins = np.zeros(grid_size, dtype=object)
            inner_occu = np.zeros(grid_size, dtype=bool)
            outer_bins = np.zeros(grid_size, dtype=object)
            outer_occu = np.zeros(grid_size, dtype=bool)
            for i in range(grid_size[0]):
                for j in range(grid_size[1]):
                    inner_bins[i, j] = [] 
                    outer_bins[i, j] = []
            ### lower left of grid
            lower = (bbox[0] - grid_offset[0], bbox[1] - grid_offset[1])

            ### now, bin inner layer
            for i in range(len(inner_coords[poly_inner])):
                ### bin inner layer 
                inner_ids = (int((inner_coords[poly_inner][i][0]-lower[0])//bin_res), 
                           (int(inner_coords[poly_inner][i][1]-lower[1])//bin_res))
                inner_bins[inner_ids[0], inner_ids[1]].append(i)
                inner_occu[inner_ids[0], inner_ids[1]] = True
            ### now, bin outer layer
            for i in range(len(outer)):
                ### bin outer layer
                outer_ids = (int((outer[i][0]-lower[0])//bin_res), (int(outer[i][1]-lower[1])//bin_res))
                outer_bins[outer_ids[0], outer_ids[1]].append(i)
                outer_occu[outer_ids[0], outer_ids[1]] = True
            
            ### compute the list of ids of occupied outer bins...
            outer_where = np.argwhere(outer_occu == True) 

            sample_pair = None
            sample_dist = float("-inf")

            ### finally ... we can do binned nearest neighbors
            ### do root2 rings... focus on points in outer perim
            ### iterate over occupied bins in outer ring to narrow comparison
            for occu_loc in outer_where:
                ### start at the center bin and iteratively look further until we find all squares...
                ### withing 2sqrt(2) + slop of the closest point we find
                ### expand out until we find all squares with UL distance of 2sqrt(2) of the closest
                ### start with this loose upper bound because this is the furthest we can possibly iterate 
                ### ...because there are only so many bins...
                root2_dist = maxbins ### needs some work... to tighten upper bound?
                root2_set = True
                ring = 0
                ring_sqs = None
                ring_offset = None
                ### while the ring is within the upper bound
                while ring < root2_dist:
                    ### get shortlist of occupied bins (inner perim) within the ring
                    temp = np.argwhere(inner_occu[max(occu_loc[0]-ring, 0): min(occu_loc[0]+ring+1, grid_size[0] - 1), 
                                            max(occu_loc[1]-ring, 0): min(occu_loc[1]+ring+1, grid_size[1] - 1)]==True)
                    ### if there are occupied bins in this ring (and we can assume we haven't found an occupied bin yet),
                    ### ... we can cap the number of rings with a loose-ish upper bound of the L-inf distance
                    ### in which we could find a point with a smaller L2-distance
                    if len(temp) > 0 and root2_set:
                        ### if this is a spot we don't need to worry about bounds
                        if spot_flag:
                            root2_dist = math.ceil(root2 * (ring + 1)) + slop
                            root2_set = False
                        ### otherwise, need to actually check this combination stays within bounds 
                        ### ... before we know we have found a valid pair
                        else:
                            for k in range(len(outer_bins[occu_loc[0], occu_loc[1]])):
                                for i in range(len(temp)):
                                    temp_box = inner_bins[temp[i][0] + max(occu_loc[0]-ring, 0), 
                                                            temp[i][1] + max(occu_loc[1]-ring, 0)]
                                    for l in range(len(temp_box)):
                                        ## s1c[bins_1[occu_loc[0], occu_loc[1]][k]]
                                        if buffer_poly.contains(LineString([(inner_coords[poly_inner][temp_box[l]]),
                                                                                    outer[outer_bins[occu_loc[0], occu_loc[1]][k]]])):
                                            root2_dist = math.ceil(root2 * (ring + 1)) + slop
                                            root2_set = False
                                            break
                                    if not root2_set:
                                        break
                                if not root2_set:
                                    break
                    ### if we have reached the upper bound and have found points...
                    ### set nearest-neighbor-check shortlist to temp
                    if ring+1 >= root2_dist and not root2_set:
                        ring_sqs = temp
                        ### need to compute lower coords of ring... because temp locations are relative to this 
                        ring_offset = [max(occu_loc[0]-ring, 0), max(occu_loc[1]-ring, 0)]
                    elif not root2_set:
                        ring = root2_dist-2
                    ring += 1
                ### for each outer box (occu_loc) we have a series of points (len(bins_1[][])) -- 
                ###     for each point we have a series of inner boxes (ring_sqs[][][i/j])
                ###         for each inner box we have a series of points --- measure distances and find min
                spair = None
                sdist = float("-inf")
                n_checks = 0
                ### iterate over every point in this outer-perim bin
                for k in range(len(outer_bins[occu_loc[0], occu_loc[1]])):
                    tpair = None
                    tdist = float("inf")
                    ### iterate over every occupied bin found in the ring process..
                    for i in range(len(ring_sqs)):
                        ### get list of points from occupied bin i
                        inner_box = inner_bins[ring_sqs[i][0] + ring_offset[0], ring_sqs[i][1] + ring_offset[1]]
                        ### iterate over every point in this box
                        for l in range(len(inner_box)):
                            ### don't need to worry about bounds if we think its a spot...
                            if spot_flag or buffer_poly.contains(LineString([outer[outer_bins[occu_loc[0], occu_loc[1]][k]],
                                                                                     inner_coords[poly_inner][inner_box[l]]])):
                                ### we can NOW compare distances between every point in the central outer-perim bin and every
                                ### point in this inner-perim bin
                                ### distfunc(s1c[bins_1[occu_loc[0], occu_loc[1]][k]], s0c[temp_box_0[l]], distparams)
                                temp_dist = compute_dist(outer[outer_bins[occu_loc[0], occu_loc[1]][k]], inner_coords[poly_inner][inner_box[l]])
                                ### check if we have found a new shortest pair...
                                ### this results in a shortest-pair combination for every point in this outer bin
                                if temp_dist < tdist:
                                    tpair = (inner_coords[poly_inner][inner_box[l]], outer[outer_bins[occu_loc[0], occu_loc[1]][k]])
                                    tdist = temp_dist
                            n_checks += 1
                    ### now having compared all inner samples to this point...
                    ### we can compare the shortest-pair distance we found to the other ones in this bin...
                    ### and hold on to the longest of them...
                    if tdist > sdist:
                        spair = tpair
                        sdist = tdist
                ### now we can compare the longest shortest-pair combo we found from this bin to the ones
                ### found in other bins...
                ### TODO -- some way to account for ties...
                if sdist > sample_dist:
                    sample_dist = sdist
                    sample_pair = spair
            ### at this point we have the furthest nearest-neighbor for this polygon pair
            if sample_dist < poly_min_dist:
                poly_min_dist = sample_dist
                poly_coordpair = sample_pair
                poly_pair = (poly_inner, poly_outer)
        ### aggregate over all outer polygons
        result_dist.append(poly_min_dist)
        result_coord_pair.append(poly_coordpair)
        result_poly_pair.append(poly_pair)
    max_loc = np.argmax(result_dist)
    maximum_distance = result_dist[max_loc]
    max_dist_origin = result_coord_pair[max_loc][0]
    max_dist_destination = result_coord_pair[max_loc][1]
    return maximum_distance, max_dist_origin, max_dist_destination, outer_coords

### fire distance computations
### direct (as crow flies) distance
def compute_dist(a, b):
    return math.sqrt(((a[0] - b[0]) ** 2) + ((a[1] - b[1]) ** 2))