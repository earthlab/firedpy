# -*- coding: utf-8 -*-
from bs4 import BeautifulSoup
from collections import OrderedDict
import datetime as dt
import ftplib
import gc
import geopandas as gpd
from glob import glob
from http.cookiejar import CookieJar
from netCDF4 import Dataset
import numpy as np
import os
from osgeo import gdal, ogr, osr
import pandas as pd
import rasterio
from rasterio.merge import merge
from shapely.geometry import Point, Polygon, MultiPolygon
import sys
import time
from tqdm import tqdm
import urllib.request as urllib2
import warnings
import yaml

# Supress depreciation warning for dask when importing xarray (annoying)
yaml.warnings({"YAMLLoadWarning": False})
warnings.filterwarnings("ignore", category=FutureWarning)
import xarray as xr
pd.options.mode.chained_assignment = None  # default="warn"


# Functions
def buildEvents(dest, data_dir, tiles, spatial_param=5, temporal_param=11):
    """
    Use the EventGrid class to classify events tile by tile and then merge
    them all together for a seamless set of wildfire events.
    """
    if not(os.path.exists(os.path.dirname(dest))):
        os.makedirs(os.path.dirname(dest))

    # Get the requested files and build an empty data frame
    files = []
    for t in tiles:
        path = os.path.join(data_dir, "rasters/burn_area/netcdfs/" + t + ".nc")
        files.append(path)
    tile_list = []
    columns = ["id", "date", "x", "y", "duration", "edge", "tile"]
    df = pd.DataFrame(columns=columns)
    base = dt.datetime(1970, 1, 1)
    files.sort()

    # Use the first file to extract some more parameters for later
    builder = EventGrid(nc_path=files[0],
                        data_dir=data_dir,
                        spatial_param=spatial_param,
                        temporal_param=temporal_param)
    geom = builder.data_set.crs.geo_transform
    res = geom[1]
    sp_buf = builder.spatial_param * res

    # Loop through each netcdf file and build individual tile events
    for file in files:
        tile_id = file[-9:-3]
        if os.path.exists(os.path.join(data_dir, "tables/events/" + tile_id +
                                       ".csv")):
            print(tile_id  + " event table exists, skipping...")
        elif not os.path.exists(
                os.path.join(data_dir, "rasters/burn_area/netcdfs/" + tile_id +
                             ".nc")):
            pass
        else:
            print("\n" + tile_id)

            # Create a new event object
            builder = EventGrid(nc_path=file, data_dir=data_dir,
                                spatial_param=5, temporal_param=11)

            # Classify event perimeters
            perims = builder.get_event_perimeters_3d()

            # Remove empty perimeters
            perims = [p for p in perims if type(p.coords[0]) is not str]
            tile_list.append(perims)

            # Extract just the event ID, days, and x,y MODIS coordinates
            plist = [(p.get_event_id(), p.coords) for p in perims]

            # Identify if it"s an edge case, so either x or y is within 5 cells
            maxys = builder.data_set["y"].data[:builder.spatial_param]
            minys = builder.data_set["y"].data[-builder.spatial_param:]
            maxxs = builder.data_set["x"].data[-builder.spatial_param:]
            minxs = builder.data_set["x"].data[:builder.spatial_param]
            yedges = list(maxys) + list(minys)
            xedges = list(maxxs) + list(minxs)

            # Create an empty data frame
            print("Building data frame...")
            events = []
            coords = []
            edges = []
            ys = []
            xs = []
            dates = []
            durations = []
            for p in plist:
                coord = [list(c) for c in p[1]]
                edge = [edgeCheck(yedges, xedges, c, sp_buf) for c in coord]
                if any(edge):
                    edge = [True for e in edge]
                event = list(np.repeat(p[0], len(coord)))
                y = [c[0] for c in coord]
                x = [c[1] for c in coord]
                date = [base + dt.timedelta(c[2]) for c in coord]
                duration = (max(date) - min(date)).days + 1
                duration = list(np.repeat(duration, len(coord)))
                events.append(event)
                coords.append(coord)
                edges.append(edge)
                ys.append(y)
                xs.append(x)
                dates.append(date)
                durations.append(duration)

            # Flatten each list of lists
            events = flttn(events)
            coords = flttn(coords)
            edges = flttn(edges)
            ys = flttn(ys)
            xs = flttn(xs)
            dates = flttn(dates)
            durations = flttn(durations)
            edf = pd.DataFrame(OrderedDict({"id": events, "date": dates, 
                                            "x": xs, "y": ys,
                                            "duration": durations,
                                            "edge": edges, "tile": tile_id}))
            if not os.path.exists(os.path.join(data_dir, "tables/events")):
                os.mkdir(os.path.join(data_dir, "tables/events"))
            edf.to_csv(
                os.path.join(data_dir, "tables/events/" + tile_id + ".csv"),
                index=False)

    # Clear memory
    gc.collect()

    # Now read in the event data frames
    print("Reading saved event tables back into memory...")
    efiles = glob(os.path.join(data_dir, "tables/events/*csv"))
    efiles = [e for e in efiles if e[-10:-4] in tiles]
    edfs = [pd.read_csv(e) for e in efiles]  # <------------------------------- Here, instead, read in list of dataframes with dask and keep as much on disk as possible?

    # Merge with existing records
    print("Concatenating event tables...")
    df = pd.concat(edfs)
    def toDays(date, base):
        date = dt.datetime.strptime(date, "%Y-%m-%d")
        delta = (date - base)
        days = delta.days
        return days

    print("Creating unique ids...")
    df["id"] = df["tile"] + "_" + df["id"].astype(str)

    print("Converting days since 1970 to dates...")
    df["days"] = df["date"].apply(toDays, base=base)

    # Cut the edge events out into a separate df
    print("Separating tile edge events from center events...")
    edges = df[df["edge"] == True]
    not_edges = df[df["edge"] == False]

    # Merge where needed
    print("Merging edge-case tile events...")
    eids = list(edges["id"].unique())
    for iden in tqdm(eids, position=0):
        # Split, one vs all
        edf = edges[edges["id"] == iden]
        edf2 = edges[edges["id"] != iden]
        days = edf["days"]

        # Sometimes these are empty?
        try:
            d1 = min(days)
            d2 = max(days)
        except:
            pass

        # If there aren"t events close enough in time the list will be empty
        edf2 = edf2[(abs(edf2["days"] - d1) < 11) | (abs(edf2["days"] - d2) < 11)]
        eids2 = list(edf2["id"].unique())

        # If there are event close in time, are they close in space?
        for iden2 in eids2:
            edf2 = edges[edges["id"] == iden2]
            ydiffs = [y - edf2["y"].values for y in edf["y"].values]
            xdiffs = [x - edf2["x"].values for x in edf["x"].values]
            ychecks = [spCheck(yds, sp_buf) for yds in ydiffs]
            xchecks = [spCheck(xds, sp_buf) for xds in xdiffs]
            checks = [ychecks[i] * xchecks[i] for i in range(len(ychecks))]
            if any(checks):
                # Merge events! Merge into the earliest event
                d12 = edf2["days"].min()
                if d1 < d12:
                    edges["id"][edges["id"] == iden2] = iden
                else:
                    edges["id"][edges["id"] == iden] = iden2

    # Concatenate edge df back into main df
    print("Recombining edge and center cases...")
    df = pd.concat([not_edges, edges])

    # Reset id values in chronological order
    print("Resetting ids in chronological order..")
    df["first"] = df.groupby("id").days.transform("min")
    firsts = df[["id", "first"]].drop_duplicates()
    firsts = firsts.sort_values("first")
    firsts["new_id"] = range(1, firsts.shape[0] + 1)
    idmap = dict(zip(firsts["id"], firsts["new_id"]))
    df["id"] = df["id"].map(idmap)
    df = df.sort_values("id")

    # Let"s add a duration and detections field - might not be necessary
    detections = df.groupby("id").size().to_frame("count").reset_index()
    detsmap = dict(zip(detections["id"], detections["count"]))
    df["detections"] = df["id"].map(detsmap)
    durations = df.groupby("id").days.apply(
                                        lambda x: max(x) - min(x)+1).to_frame()
    durmap = dict(zip(durations.index, durations["days"]))
    df["duration"] = df["id"].map(durmap)

    # put these in order
    df = df[["id", "tile", "date", "x", "y", "duration", "detections"]]

    # Finally save
    print("Saving merged event table to " + dest)
    df.to_csv(dest, index=False)


def buildPolygons(src, daily_shp_path, event_shp_path, data_dir):
    # Start the timer (seconds)
    start = time.perf_counter()

    # Make sure we have the target folders
    if not(os.path.exists(os.path.dirname(event_shp_path))):
        os.makedirs(os.path.dirname(event_shp_path))

    # Read in the event table and a reference nc file for geometric information
    print("Reading classified fire event table...")
    df = pd.read_csv(src)

    # Go ahead and create daily id (did) for later
    df["did"] = df["id"].astype(str) + "-" + df["date"]

    # Use a sample for geometries
    sample_path = glob(os.path.join(data_dir,
                                    "rasters/burn_area/netcdfs/*.nc"))[0]
    sample = Dataset(sample_path)
    crs = sample.variables["crs"]
    geom = crs.geo_transform
    proj4 = crs.proj4
    res = [geom[1], geom[-1]]

    # Filter columns, center pixel coordinates, and remove repeating pixels
    df = df[["id", "did", "date", "x", "y"]]
    df["x"] = df["x"] + (res[0]/2)
    df["y"] = df["y"] + (res[1]/2)

    # Each entry in the df gets a point object from the x and y coordinates.
    print("Converting data frame to spatial object...")
    df["geometry"] = df[["x", "y"]].apply(lambda x: Point(tuple(x)), axis=1)
    gdf = df[["id", "did", "date", "geometry"]]
    gdf = gpd.GeoDataFrame(gdf, crs=proj4, geometry=gdf["geometry"])

    # Create a circle buffer
    print("Creating buffer...")
    geometry = gdf.buffer(1 + (res[0]/2))
    gdf["geometry"] = geometry

    # Then create a square envelope around the circle
    gdf["geometry"] = gdf.envelope

    # Now add the first date of each event and merge daily event detections
    print("Dissolving polygons...")
    gdf["start_date"] = gdf.groupby("id")["date"].transform("min")  
    gdfd = gdf.dissolve(by="did", as_index=False)
    gdfd["year"] = gdfd["start_date"].apply(lambda x: x[:4])
    gdfd["month"] = gdfd["start_date"].apply(lambda x: x[5:7])
    #gdfd = gdfd[gdfd["month"].isin(["06", "07", "08", "09"])]

    # Save the daily before dissolving into event level
    print("Saving daily file to " + daily_shp_path + "...")
    gdfd.to_file(daily_shp_path, driver="GPKG")

    # Now merge into event level polygons
    gdf = gdf[["id", "start_date", "geometry"]]
    gdf = gdf.dissolve(by="id", as_index=False)

    # For each geometry, if it is a single polygon, cast as a multipolygon
    print("Converting polygons to multipolygons...")
    def asMultiPolygon(polygon):
        if type(polygon) == Polygon:
            polygon = MultiPolygon([polygon])
        return polygon
    gdf["geometry"] = gdf["geometry"].apply(asMultiPolygon)

    # Calculate perimeter length
    print("Calculating perimeter lengths...")
    gdf["final_perimeter"] = gdf["geometry"].length  # <----------------------- Check accuracy of this, QGIS is slightly different (also check adams)

    # Now save as a geopackage  # <-------------------------------------------- Should we also make a shapefile for ESRI users? Make it optional.
    print("Saving event-level file to " + event_shp_path + "...")
    gdf.to_file(event_shp_path, driver="GPKG")

    # Print the time it took
    end = time.perf_counter()
    seconds = end - start
    minutes = seconds/60
    print("Job completed in {} minutes".format(round(minutes, 2)))


def convertDates(array, year):
    """
    Convert everyday in an array to days since Jan 1 1970
    """
    def convertDate(julien_day, year):
        base = dt.datetime(1970, 1, 1)
        date = dt.datetime(year, 1, 1) + dt.timedelta(int(julien_day))
        days = date - base
        return days.days

    # Loop through each position with data and convert
    locs = np.where(array > 0)
    ys = locs[0]
    xs = locs[1]
    locs = [[ys[i], xs[i]] for i in range(len(xs))]
    for l in locs:
        y = l[0]
        x = l[1]
        array[y, x] = convertDate(array[y, x], year)

    return array


def dateRange(perimeter):
    """
    Converts days in a perimeter object since Jan 1 1970 to date strings
    """
    if len(perimeter.coords) > 0:
        base = dt.datetime(1970, 1, 1)
        days = [p[2] for p in perimeter.coords]
        day1 = (base + dt.timedelta(days=int(min(days)))).strftime("%Y-%m-%d")
    else:
        day1 = "N/A"
    return day1


def edgeCheck(yedges, xedges, coord, sp_buffer):
    """
    Identify edge cases to make merging events quicker later
    """
    y = coord[0]
    x = coord[1]
    if y in yedges:
        edge = True
    elif x in xedges:
        edge = True
    else:
        edge = False
    return edge


def flttn(lst):
    """
    Just a quick way to flatten lists of lists
    """
    lst = [l for sl in lst for l in sl]
    return lst


def maxGrowthDate(x):
    dates = x["date"].to_numpy()
    pixels = x["pixels"].to_numpy()
    loc = np.where(pixels == np.max(pixels))[0]
    d = np.unique(dates[loc])
    if len(d) > 1:
        d = ", ".join(d)
    else:
        d = d[0]
    return d


def toAcres(p, res):
    return (p*res**2) * 0.000247105


def toHa(p, res):
    return (p*res**2) * 0.0001


def toKms(p, res):
    return (p*res**2)/1000000


def mergeChecker(new_coords, full_list, temporal_param, radius):
    """
    This uses a radius for the spatial window as opposed to a square and is not
    currently being used to merge events.
    """
    t1 = np.min([c[2] for c in new_coords]) - temporal_param
    t2 = np.max([c[2] for c in new_coords]) + temporal_param
    for i in range(len(full_list)):
        old_event = full_list[i]
        old_coords = old_event[1]
        old_times = [c[2] for c in old_coords]
        time_checks = [t for t in old_times if t >= t1 and t <= t2]

        if len(time_checks) > 0:
            for coord in new_coords:
                # Check if the time coordinate is within an old event
                radii = []
                new_y = coord[0]
                new_x = coord[1]
                for oc in old_coords:
                    old_y = oc[0]
                    old_x = oc[1]
                    dy = abs(old_y - new_y)
                    dx = abs(old_x - new_x)
                    r = np.sqrt((dy ** 2) + (dx ** 2))
                    radii.append(r)
                check = [r for r in radii if r <= radius]
                if any(check):
                    return i, True
                else:
                    return i, False
            else:
                return i, False


def mode(lst):
    if len(np.unique(lst)) > 1:
        grouped_lst = [list(lst[lst == s]) for s in lst]
        counts = {len(a): a for a in grouped_lst}  # overwrites matches
        max_count = np.max(list(counts.keys()))
        mode = counts[max_count][0]
    else:
        mode = lst.unique()[0]
    return mode


def pquery(p, lc, lc_array):
    """
    Find the landcover code for a particular point (p).
    """
    row, col = lc.index(p.x, p.y)
    lc_value = lc_array[row, col]
    return lc_value


def rasterize(src, dst, attribute, resolution, crs, extent, all_touch=False,
              na=-9999):

    # Open shapefile, retrieve the layer
    src_data = ogr.Open(src)
    layer = src_data.GetLayer()

    # Use transform to derive coordinates and dimensions
    xmin = extent[0]
    ymin = extent[1]
    xmax = extent[2]
    ymax = extent[3]

    # Create the target raster layer
    cols = int((xmax - xmin)/resolution)
    rows = int((ymax - ymin)/resolution) + 1
    trgt = gdal.GetDriverByName("GTiff").Create(dst, cols, rows, 1,
                                gdal.GDT_Float32)
    trgt.SetGeoTransform((xmin, resolution, 0, ymax, 0, -resolution))

    # Add crs
    refs = osr.SpatialReference()
    refs.ImportFromWkt(crs)
    trgt.SetProjection(refs.ExportToWkt())

    # Set no value
    band = trgt.GetRasterBand(1)
    band.SetNoDataValue(na)

    # Set options
    if all_touch is True:
        ops = ["-at", "ATTRIBUTE=" + attribute]
    else:
        ops = ["ATTRIBUTE=" + attribute]

    # Finally rasterize
    gdal.RasterizeLayer(trgt, [1], layer, options=ops)

    # Close target an source rasters
    del trgt
    del src_data


def spCheck(diffs, sp_buf):
    """
    Quick function to check if events land within the spatial window.
    """
    checks = [e for e in diffs if abs(e) < sp_buf]
    if any(checks):
        check = True
    else:
        check = False
    return check


def toDays(date, base):
    """
    Convert dates to days since a base date
    """
    if type(date) is str:
        date = dt.datetime.strptime(date, "%Y-%m-%d")
        delta = (date - base)
        days = delta.days
    return days


def toRaster(array, trgt, geometry, proj, navalue=-9999):
    """
    Writes a single array to a raster with coordinate system and
    geometric information.

    trgt = target path
    proj = spatial reference system
    geom = geographic transformation
    """
    ypixels = array.shape[0]
    xpixels = array.shape[1]
    trgt = trgt.encode("utf-8")
    image = gdal.GetDriverByName("GTiff").Create(trgt, xpixels, ypixels, 1,
                                                 gdal.GDT_Float32)
    image.SetGeoTransform(geometry)
    image.SetProjection(proj)
    image.GetRasterBand(1).WriteArray(array)
    image.GetRasterBand(1).SetNoDataValue(navalue)


# Classes
class DataGetter:
    """
    Things to do/remember:
        - authorization issues, add user inputs
        - use shapefile to determine appropriate tiles
        - create a pool for parallel downloads
    """
    def __init__(self, data_dir):
        self.data_dir = data_dir
        self.date = dt.datetime.today().strftime("%m-%d-%Y")
        self.createPaths()
        self.cpus = os.cpu_count()
        self.modis_template_path = os.path.join(data_dir, "rasters/")
        self.modis_template_file_root = "mosaic_template.tif"
        self.landcover_path = os.path.join(data_dir, "rasters/landcover")
        self.landcover_mosaic_path = os.path.join(data_dir,
                                                  "rasters/landcover/mosaics")
        self.landcover_file_root = "us_lc_mosaic_"
        self.nc_path = os.path.join(data_dir, "rasters/burn_area/netcdfs")
        self.hdf_path = os.path.join(data_dir, "rasters/burn_area/hdfs")
        self.tiles = ["h08v04", "h09v04", "h10v04", "h11v04", "h12v04",
                      "h13v04", "h08v05", "h09v05", "h10v05", "h11v05",
                      "h12v05", "h08v06", "h09v06", "h10v06", "h11v06"]
        print("Project Folder: " + data_dir)

    def createPaths(self):
        sub_folders = ["rasters", "shapefiles", "shapefiles/ecoregion",
                       "tables", "rasters/burn_area", "rasters/burn_area/hdfs",
                       "rasters/landcover", "rasters/landcover/mosaics/",
                       "rasters/ecoregion"]
        folders = [os.path.join(self.data_dir, sf) for sf in sub_folders]
        for f in folders:
            if not os.path.exists(f):
                os.mkdir(f)

    def downloadBA(self, file, ftp):
        missing = []
        tile = file[17:23]
        folder = os.path.join(self.hdf_path, tile)
        if not os.path.exists(folder):
            os.mkdir(folder)
        trgt = os.path.join(folder, file)
        if not os.path.exists(trgt):
            try:
                with open(trgt, "wb") as dst:
                    ftp.retrbinary("RETR %s" % file, dst.write, 102400)
            except ftplib.all_errors as e:
                missing.append(file)
                print("FTP Transfer Error: ", e)
        return missing

    def getBurns(self):
        """
        This will download the MODIS burn event data set tiles and create a
        singular mosaic to use as a template file for coordinate reference
        information and geometries.

        User manual:
            http://modis-fire.umd.edu/files/MODIS_C6_BA_User_Guide_1.2.pdf

        FTP:
            ftp://fire:burnt@fuoco.geog.umd.edu/gfed4/MCD64A1/C6/
        """
        # Check in to the site
        ftp = ftplib.FTP("fuoco.geog.umd.edu", user="fire", passwd="burnt")
        ftp.cwd("MCD64A1/C6")

        # Use specified tiles or...download all tiles if the list is empty
        if self.tiles[0].lower() != "all":
            tiles = self.tiles
        else:
            tiles = ftp.nlst()
            tiles = [t for t in tiles if "h" in t]       

        # Download the available files and catch failed downloads
        missings = []
        for tile in tiles:
            if not os.path.exists(
                     os.path.join(self.data_dir, "rasters/burn_area/netcdfs/" + 
                                  tile + ".nc")):
                print("Downloading/Checking hdf files for " + tile)
                ftp_folder =  "/MCD64A1/C6/" + tile
                try:
                    ftp.cwd(ftp_folder)
                    hdfs = ftp.nlst()
                    hdfs = [h for h in hdfs if ".hdf" in h]
                    for h in tqdm(hdfs, position=0):
                        trgt = os.path.join(self.hdf_path, tile, h)
    
                        # Check if the file already exists and works
                        if os.path.exists(trgt):
                            missing = []
                        else:
                            missing = self.downloadBA(h, ftp)
    
                        # If the download failed, add to missing
                        if len(missing) > 0:
                            missings = missings + missing
    
                        # Even if it didn"t fail, make sure it works
                        try:
                            gdal.Open(trgt).GetSubDatasets()[0][0]
                        except:
                            print("Bad file, removing and trying again...")
                            missings = missings + [h]
                            os.remove(trgt)
                except:
                    pass

        # Now try again for the missed files
        if len(missings) > 0:
            print("Missed Files: \n" + str(missings))
            print("trying again...")
            for file in missings:
                tile = file[17:23]
                ftp_folder =  "/MCD64A1/C6/" + tile
                ftp.cwd(ftp_folder)
                trgt = os.path.join(self.hdf_path, tile, h)
                try:
                    with open(trgt, "wb") as dst:
                        ftp.retrbinary("RETR %s" % file, dst.write, 102400)
                except ftplib.all_errors as e:
                    missing.append(file)
                    print("FTP Transfer Error: ", e)

                try:
                    gdal.Open(trgt).GetSubDatasets()[0][0]
                    missings.remove(file)
                except Exception as e:
                    print(e)

        # If that doesn"t get them all, give up and let the user handle it.
        if len(missings) > 0:
            print("There are still " + str(len(missings)) + " missed files.")
            print("Try downloading these files manually: ")
            for m in missings:
                print(m)

        # And quit the ftp service
        ftp.quit()

        # Build the netcdfs here
        tile_files = {}
        for tid in tiles:
            files = glob(os.path.join(self.hdf_path, tid, "*hdf"))
            tile_files[tid] = files

        # Merge one year into a reference mosaic
        if not os.path.exists(self.modis_template_path):
            folders = glob(os.path.join(self.hdf_path, "*"))
            file_groups = [glob(os.path.join(f, "*hdf")) for f in folders]
            for f in file_groups:
                f.sort()
            files = [f[0] for f in file_groups]
            dss = [rasterio.open(f).subdatasets[0] for f in files]
            tiles = [rasterio.open(d) for d in dss]
            mosaic, transform = merge(tiles)
            crs = tiles[0].meta.copy()
            template_path = os.path.join(self.modis_template_path,
                                         self.modis_template_file_root)
            crs.update({"driver": "GTIFF",
                        "height": mosaic.shape[1],
                        "width": mosaic.shape[2],
                        "transform": transform})
            with rasterio.open(template_path, "w", **crs) as dest:
                dest.write(mosaic)

        # Build one netcdf per tile
        for tid in tiles:
            files = tile_files[tid]
            if len(files) > 0:
                try:
                    self.buildNCs(files)
                except Exception as e:
                    file_name = os.path.join(self.nc_path, tid + ".nc")
                    print("Error on tile " + tid + ": " + str(e))
                    print("Removing " + file_name + " and moving on.")
                    os.remove(file_name)

    def getLandcover(self):
        """
        A method to download and process landcover data from NASA"s Land
        Processes Distributed Active Archive Center, which is an Earthdata
        thing. You"ll need register for a username and password, but that"s
        free. Fortunately, there is a tutorial on how to get this data:
            
        https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+
        Python
        
        sample citation for later:
            ASTER Mount Gariwang image from 2018 was retrieved from
            https://lpdaac.usgs.gov, maintained by the NASA EOSDIS Land
            Processes Distributed Active Archive Center (LP DAAC) at the USGS
            Earth Resources Observation and Science (EROS) Center, Sioux Falls,
            South Dakota. 2018, https://lpdaac.usgs.gov/resources/data-action/
            aster-ultimate-2018-winter-olympics-observer/.
        """
        # Use specified tiles or...
        if self.tiles[0].lower() != "all":
            tiles = self.tiles

        # ...download all tiles if the list is empty
        else:
            # Check in to the burn data site
            ftp = ftplib.FTP("fuoco.geog.umd.edu")
            ftp.login("fire", "burnt")
            ftp.cwd("/MCD64A1/C6/")
            tiles = ftp.nlst()
            tiles = [t for t in tiles if "h" in t]

        # We need years in strings  <------------------------------------------ Parameterize this
        years = [str(y) for y in range(2001, 2017)]

        # Access
        print("Retrieving land cover rasters from NASA's Earthdata service...")
        print("Register at the link below to obtain a username and password:")
        print("https://urs.earthdata.nasa.gov/")
        username = input("Enter NASA Earthdata User Name: ")
        password = input("Enter NASA Earthdata Password: ")
        pw_manager = urllib2.HTTPPasswordMgrWithDefaultRealm()
        pw_manager.add_password(None, "https://urs.earthdata.nasa.gov",
                                username, password)
        cookiejar = CookieJar()
        opener = urllib2.build_opener(urllib2.HTTPBasicAuthHandler(pw_manager),
                                      urllib2.HTTPCookieProcessor(cookiejar))
        urllib2.install_opener(opener)

        # Land cover data from earthdata.nasa.gov
        for year in years:
            print("Retrieving landcover data for " + year)
            url = ("https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/" + year +
                   ".01.01/")
            r = urllib2.urlopen(url)
            soup = BeautifulSoup(r, features="lxml", 
                                 from_encoding=r.info().get_param("charset")
                                 )
            names = [link["href"] for link in soup.find_all("a", href=True)]
            names = [f for f in names if "hdf" in f and f[17:23] in tiles]
            links = [url + l for l in names]
            for i in tqdm(range(len(links)), position=0):
                if not os.path.exists(os.path.join(self.landcover_path, year)):
                    os.mkdir(os.path.join(self.landcover_path, year))
                path = os.path.join(self.landcover_path, year, names[i])
                if not os.path.exists(path):
                    request = urllib2.Request(links[i])
                    with open(path, "wb") as file:
                        response = urllib2.urlopen(request).read()
                        file.write(response)

        # Now process these tiles into yearly geotiffs.
        if not os.path.exists(self.landcover_mosaic_path):
            os.mkdir(self.landcover_mosaic_path)
            os.mkdir(os.path.join(self.landcover_mosaic_path, "wgs"))
        for year in years:
            print("Stitching together landcover tiles for year " + year)
            lc_tiles = glob(os.path.join(self.landcover_path, year, "*hdf"))
            dss = [rasterio.open(f).subdatasets[0] for f in lc_tiles]           
            tiles = [rasterio.open(d) for d in dss]
            mosaic, transform = merge(tiles)
            crs = tiles[0].meta.copy()
            crs.update({"driver": "GTIFF",
                        "height": mosaic.shape[1],
                        "width": mosaic.shape[2],
                        "transform": transform})
            file = self.landcover_file_root + year + ".tif"
            path = os.path.join(self.landcover_mosaic_path, file)
            with rasterio.open(path, "w", **crs) as dest:
                dest.write(mosaic)


    def getShapes(self):
        """
        Just to grab some basic shapefiles needed for calculating statistics.
        """
        if not os.path.exists(os.path.join(self.data_dir, "shapefiles")):
            os.mkdir(os.path.join(self.data_dir, "shapefiles"))

        # Variables
        conus_states = ["WV", "FL", "IL", "MN", "MD", "RI", "ID", "NH", "NC",
                        "VT", "CT", "DE", "NM", "CA", "NJ", "WI", "OR", "NE",
                        "PA", "WA", "LA", "GA", "AL", "UT", "OH", "TX", "CO",
                        "SC", "OK", "TN", "WY", "ND", "KY", "VI", "ME", "NY",
                        "NV", "MI", "AR", "MS", "MO", "MT", "KS", "IN", "SD",
                        "MA", "VA", "DC", "IA", "AZ"]
        modis_crs = ("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 " +
                     "+b=6371007.181 +units=m +no_defs")

        # MODIS Sinusoial World Grid
        if not os.path.exists(
                os.path.join(self.data_dir,
                             "shapefiles/modis_world_grid.shp")):
            print("Downloading MODIS Sinusoidal Projection Grid...")
            src = ("http://book.ecosens.org/wp-content/uploads/2016/06/" +
                   "modis_grid.zip")
            modis = gpd.read_file(src)
            modis.crs = modis_crs
            modis.to_file(os.path.join(self.data_dir,
                                       "shapefiles/modis_world_grid.shp"))

        # Contiguous United States - WGS84
        if not os.path.exists(os.path.join(self.data_dir,
                                           "shapefiles/conus.shp")):
            print("Downloading US state shapefile from the Census Bureau...")
            usa = gpd.read_file("http://www2.census.gov/geo/tiger/GENZ2016/" +
                                "shp/cb_2016_us_state_20m.zip")
            conus = usa[usa["STUSPS"].isin(conus_states)]
            conus.crs = {"init": "epsg:4326", "no_defs": True}
            conus.to_file(os.path.join(self.data_dir, "shapefiles/conus.shp"))

        # Contiguous United States - MODIS Sinusoidal
        if not os.path.exists(os.path.join(self.data_dir, 
                                           "shapefiles/conus_modis.shp")):
            print("Reprojecting state shapefile to MODIS Sinusoidal...")
            conus = gpd.read_file(os.path.join(self.data_dir,
                                               "shapefiles/conus.shp"))
            modis_conus = conus.to_crs(modis_crs)
            modis_conus.to_file(os.path.join(self.data_dir, 
                                             "shapefiles/conus_modis.shp"))

        # Level III Omernick Ecoregions - USGS North American Albers
        if not os.path.exists(
                os.path.join(self.data_dir,
                             "shapefiles/ecoregion/us_eco_l3.shp")):
            print("Downloading Omernick Level III Ecoregions from the USGS...")
            eco_l3 = gpd.read_file("ftp://ftp.epa.gov/wed/ecoregions/us/" +
                                   "us_eco_l3.zip")
            eco_l3.crs = {"init": "epsg:5070"}
            eco_l3.to_file(os.path.join(self.data_dir, 
                                        "shapefiles/ecoregion/us_eco_l3.shp"))
            eco_l3 = eco_l3.to_crs(modis_crs)
            eco_l3.to_file(
                    os.path.join(self.data_dir,
                                 "shapefiles/ecoregion/us_eco_l3_modis.shp"))
            eco_ref = eco_l3[["US_L3CODE", "NA_L3NAME", "NA_L2NAME",
                              "NA_L1NAME"]].drop_duplicates()
            def cap(string):
                strings = string.split()
                strings = [s.lower() if s != "USA" else s for s in strings]
                caps = [s.capitalize()  if s != "USA" else s for s in strings]
                return " ".join(caps)
            eco_ref["NA_L2NAME"] = eco_ref["NA_L2NAME"].apply(cap)
            eco_ref["NA_L1NAME"] = eco_ref["NA_L1NAME"].apply(cap)
            eco_ref.to_csv(os.path.join(self.data_dir, "tables/eco_refs.csv"),
                           index=False)

        # Rasterize Level Omernick Ecoregions - WGS 84
        if not os.path.exists(
                os.path.join(self.data_dir,
                             "rasters/ecoregion/us_eco_l3_modis.tif")):

            # We need something with the correct geometry
            src = os.path.join(self.data_dir,
                               "shapefiles/ecoregion/us_eco_l3_modis.shp")
            dst = os.path.join(self.data_dir,
                               "rasters/ecoregion/us_eco_l3_modis.tif")
            extent_template_file = os.path.join(
                    self.data_dir, "shapefiles/modis_world_grid.shp")

            # Getting the extent regardless of existing files from other runs
            template1 = gpd.read_file(extent_template_file)            
            template1["h"] = template1["h"].apply(lambda x: "{:02d}".format(x))
            template1["v"] = template1["v"].apply(lambda x: "{:02d}".format(x))
            template1["tile"] = "h" + template1["h"] + "v" +  template1["v"]
            template1 = template1[template1["tile"].isin(self.tiles)]

            # We can use this to query which tiles are needed for coordinates
            bounds = template1.geometry.bounds
            minx = min(bounds["minx"])
            miny = min(bounds["miny"])
            maxx = max(bounds["maxx"])
            maxy = max(bounds["maxy"])
            minx_tile = template1["tile"][bounds["minx"] == minx].iloc[0]
            miny_tile = template1["tile"][bounds["miny"] == miny].iloc[0]
            maxx_tile = template1["tile"][bounds["maxx"] == maxx].iloc[0]
            maxy_tile = template1["tile"][bounds["maxy"] == maxy].iloc[0]
            extent_tiles = [minx_tile, miny_tile, maxx_tile, maxy_tile]

            # If these aren"t present, I say just go ahead and download
            exts = []
            for et in extent_tiles:
                folder = os.path.join(self.data_dir, "rasters/burn_area/hdfs",
                                      et)
                if not os.path.exists(folder):
                    self.getBurns()
                file = glob(os.path.join(folder, "*hdf"))[0]
                file_pointer = gdal.Open(file)
                dataset_pointer = file_pointer.GetSubDatasets()[0][0]
                ds = gdal.Open(dataset_pointer)
                geom = ds.GetGeoTransform()
                ulx, xres, xskew, uly, yskew, yres = geom
                lrx = ulx + (ds.RasterXSize * xres)
                lry = uly + (ds.RasterYSize * yres) + yres  # <---------------- Off by one cell
                exts.append([ulx, lry, lrx, uly])

            extent = [exts[0][0], exts[1][1], exts[2][2], exts[3][3]]
            wkt = ds.GetProjection()
            attribute = "US_L3CODE"
            rasterize(src, dst, attribute, xres, wkt, extent)

    def shapeToTiles(self, shp_path):
        """
        Set or reset the tile list using a shapefile. Where shapes intersect
        with the modis sinusoidal grid determines which tiles to use.
        """
        source = gpd.read_file(shp_path)

        modis_crs = {'proj': 'sinu', 'lon_0': 0, 'x_0': 0, 'y_0': 0,
                     'a': 6371007.181, 'b': 6371007.181, 'units': 'm',
                     'no_defs': True}

        # Attempt to project to modis sinusoidal
        try:
            source = source.to_crs(modis_crs)
        except Exception as e:
            print("Error: " + str(e))
            print("Failed to reproject file, ensure a coordinate reference " +
                  "system is specified.")

        # Attempt to read in the modis grid and download it if not available
        try:
            modis_grid = gpd.read_file(
                    os.path.join(self.data_dir,
                                 "shapefiles/modis_world_grid.shp"))
        except:
            modis_grid = gpd.read_file("http://book.ecosens.org/wp-content/" +
                                       "uploads/2016/06/modis_grid.zip")
            modis_grid.to_file(os.path.join(self.data_dir,
                                            "shapefiles/modis_world_grid.shp"))

        # Left join shapefiles with source shape as the left
        shared = gpd.sjoin(source, modis_grid, how="left").dropna()
        shared["h"] = shared["h"].apply(lambda x: "h{:02d}".format(int(x)))
        shared["v"] = shared["v"].apply(lambda x: "v{:02d}".format(int(x)))
        shared["tile"] = shared["h"] + shared["v"]
        tiles = pd.unique(shared["tile"].values)
        self.tiles = tiles


    def buildNCs(self, files):
        """
        Take in a time series of files for the MODIS burn detection dataset and
        create a singular netcdf file.
        """
        savepath = self.nc_path

        # Check that the target folder exists, agian.
        if not os.path.exists(savepath):
            os.mkdir(savepath)

        # Set file names
        files.sort()
        names = [os.path.split(f)[-1] for f in files]
        names = [f.split(".")[2] + "_" + f.split(".")[1][1:] for f in names]
        tile_id = names[0].split("_")[0]
        file_name = os.path.join(savepath, tile_id + ".nc")

        # Skip if it exists already
        if os.path.exists(file_name):
            print(tile_id + " netCDF file exists, skipping...")
        else:
            # Use a sample to get geography information and geometries
            print("Building netcdf for tile " + tile_id)
            sample = files[0]
            ds = gdal.Open(sample).GetSubDatasets()[0][0]
            hdf = gdal.Open(ds)
            geom = hdf.GetGeoTransform()
            proj = hdf.GetProjection()
            data = hdf.GetRasterBand(1)
            crs = osr.SpatialReference()
    
            # Get the proj4 string usign the WKT
            crs.ImportFromWkt(proj)
            proj4 = crs.ExportToProj4()
    
            # Use one tif (one array) for spatial attributes
            array = data.ReadAsArray()
            ny, nx = array.shape
            xs = np.arange(nx) * geom[1] + geom[0]
            ys = np.arange(ny) * geom[5] + geom[3]

            # Todays date for attributes
            todays_date = dt.datetime.today()
            today = np.datetime64(todays_date)
    
            # Create Dataset
            nco = Dataset(file_name, mode="w", format="NETCDF4", clobber=True)
    
            # Dimensions
            nco.createDimension("y", ny)
            nco.createDimension("x", nx)
            nco.createDimension("time", None)
    
            # Variables
            y = nco.createVariable("y",  np.float64, ("y",))
            x = nco.createVariable("x",  np.float64, ("x",))
            times = nco.createVariable("time", np.int64, ("time",))
            variable = nco.createVariable("value",np.int16,
                                          ("time", "y", "x"),
                                          fill_value=-9999, zlib=True)
            variable.standard_name = "day"
            variable.long_name = "Burn Days"
    
            # Appending the CRS information 
            # Check "https://cf-trac.llnl.gov/trac/ticket/77"
            crs = nco.createVariable("crs", "c")
            variable.setncattr("grid_mapping", "crs")
            crs.spatial_ref = proj4
            crs.proj4 = proj4
            crs.geo_transform = geom
            crs.grid_mapping_name = "sinusoidal"
            crs.false_easting = 0.0
            crs.false_northing = 0.0
            crs.longitude_of_central_meridian = 0.0
            crs.longitude_of_prime_meridian = 0.0
            crs.semi_major_axis = 6371007.181
            crs.inverse_flattening = 0.0

            # Coordinate attributes
            x.standard_name = "projection_x_coordinate"
            x.long_name = "x coordinate of projection"
            x.units = "m"
            y.standard_name = "projection_y_coordinate"
            y.long_name = "y coordinate of projection"
            y.units = "m"
    
            # Other attributes
            nco.title = "Burn Days"
            nco.subtitle = "Burn Days Detection by MODIS since 1970."
            nco.description = "The day that a fire is detected."
            nco.date = pd.to_datetime(str(today)).strftime("%Y-%m-%d")
            nco.projection = "MODIS Sinusoidal"
            nco.Conventions = "CF-1.6"
    
            # Variable Attrs
            times.units = "days since 1970-01-01"
            times.standard_name = "time"
            times.calendar = "gregorian"
            datestrings = [f[-7:] for f in names]
            dates = []
            for d in datestrings:
                year = dt.datetime(year=int(d[:4]), month=1, day=1)
                date = year + dt.timedelta(int(d[4:]))
                dates.append(date)
            deltas = [d - dt.datetime(1970, 1, 1) for d in dates]
            days = np.array([d.days for d in deltas])
    
            # Write dimension data
            x[:] = xs
            y[:] = ys
            times[:] = days
    
            # One file a time, write the arrays
            tidx = 0
            for f in tqdm(files, position=0, file=sys.stdout):
                ds = gdal.Open(f).GetSubDatasets()[0][0]
                hdf = gdal.Open(ds)
                data = hdf.GetRasterBand(1)
                array = data.ReadAsArray()
                year = int(f[-36: -32])
                array = convertDates(array, year)
                try:
                    variable[tidx, :, :] = array
                except:
                    print(f + ": failed, probably had wrong dimensions, " +
                          "inserting a blank array in its place.")
                    blank = np.zeros((ny, nx))
                    variable[tidx, :, :] = blank
                tidx += 1
    
            # Done
            nco.close()

    def buildTiffs(self):
        """
        This will take the original HDF4 files retrieved from get Burns,
        convert the data (in julien days) to days since 1970 and save them as
        individual tiffs in "data/rasters/burn_area/geotiffs"
        """
        # Let"s start with one tile
        tiles = self.tiles
        tile_files = glob(os.path.join(self.hdf_path, "*"))
        tile_files.sort()

        # Loop through everything, convert dates, and save to tiff
        for t in tile_files:
            files = glob(os.path.join(t, "*hdf"))
            files.sort()
            for f in tqdm(files, position=0):
                if not os.path.exists(f):
                    # Get data and attributes
                    year = int(f[44:48])
                    day = f[48:51]
                    ds = gdal.Open(f).GetSubDatasets()[0][0]
                    hdf = gdal.Open(ds)
                    geometry = hdf.GetGeoTransform()
                    proj = hdf.GetProjection()
                    data = hdf.GetRasterBand(1)
                    array = data.ReadAsArray()

                    # Convert dates in array
                    array = convertDates(array, year)

                    # Save as a geotiff...see if we get the same output
                    tile = os.path.basename(t)
                    tfolder = os.path.join(self.data_dir, 
                                           "rasters/burn_area/geotiffs/",
                                           tile)
                    if not os.path.exists(tfolder):
                        os.mkdir(tfolder)
                    tfile = tile + "_" + str(year) + day + ".tif"
                    trgt = os.path.join(tfolder, tfile)
                    toRaster(array, trgt, geometry, proj)

        # Merge one year into a reference mosaic
        if not os.path.exists(self.modis_template_path):
            folder = os.path.join(self.data_dir, "rasters/burn_area/geotiffs")
            folders = glob(os.path.join(folder, "*"))
            file_groups = [glob(os.path.join( f, "*tif")) for f in folders]
            for f in file_groups:
                f.sort()
            files = [f[0] for f in file_groups]
            tiles = [rasterio.open(f) for f in files]
            mosaic, transform = merge(tiles)
            crs = tiles[0].meta.copy()
            crs.update({"driver": "GTIFF",
                        "height": mosaic.shape[1],
                        "width": mosaic.shape[2],
                        "transform": transform})
            with rasterio.open(self.modis_template_path, "w", **crs) as dest:
                dest.write(mosaic)      


class EventPerimeter:
    def __init__(self, event_id, coord_list=[]):
        self.event_id = event_id
        self.merge_id = np.nan
        self.coords = []
        self.coords = self.add_coordinates(coord_list)

    def add_coordinates(self,coord_list):
        for coord in coord_list:
            self.coords.append(coord)
        return self.coords

    def get_event_id(self):
        return self.event_id

    def get_merge_id(self):
        return self.merge_id

    def get_coords(self):
        return self.coords


class EventGrid:
    """
    For a single file, however large, find sites with any burn detections
    in the study period, loop through these sites and group by the space-time
    window, save grouped events to a data frame on disk.
    """
    def __init__(self, data_dir, nc_path=("rasters/burn_area/netcdfs/"),
                 spatial_param=5, temporal_param=11, area_unit="Unknown",
                 time_unit="days since 1970-01-01"):
        self.data_dir = data_dir
        self.nc_path = os.path.join(data_dir, nc_path)
        self.spatial_param = spatial_param
        self.temporal_param = temporal_param
        self.area_unit = area_unit
        self.time_unit = time_unit
        self.event_grid = {}
        self.next_event_id = 1
        self.input_array = self.get_input_xarray()

    def get_input_xarray(self):
        burns = xr.open_dataset(self.nc_path)
        self.data_set = burns
        self.coordinates = burns.coords
        input_array = burns.value

        return input_array

    def add_event_grid(self, event_id, new_pts):
        for p in new_pts:
            entry = {p : event_id}
            self.event_grid.update(entry)

    def merge_perimeters(self, perimeters, event_id, obsolete_id):
        # set the merge id in the obsolete id
        perimeters[obsolete_id-1].merge_id = event_id
        new_pts = perimeters[obsolete_id-1].coords

        # update the event_grid and add points to event_id perimeter
        for p in new_pts:
            self.event_grid[p] = event_id
        perimeters[event_id-1].add_coordinates(new_pts)

        # set old perimeter to null
        merge_notice = "Merged with event {}".format(event_id)
        perimeters[obsolete_id-1].coords = [merge_notice, new_pts]

        return perimeters

    def get_spatial_window(self, y, x, array_dims):
        """
        Pull in the spatial window around a detected event and determine its
        shape and the position of the original point within it. Finding this
        origin point is related to another time saving step in the event
        classification procedure.
        """
        top = max(0, y - self.spatial_param)
        bottom = min(array_dims[0], y + self.spatial_param)
        left = max(0, x - self.spatial_param)
        right = min(array_dims[1], x + self.spatial_param)

        #  Derive full xarray coordinates from just the window for speed
        ydim = array_dims[0]
        xdim = array_dims[1]

        # Expand the spatial dimension
        tps = [i for i in range(self.spatial_param)]

        # There are four edge cases
        x_edges_1 = [0 + t for t in tps]
        x_edges_2 = [xdim - t for t in tps]
        y_edges_1 = [0 + t for t in tps]
        y_edges_2 = [ydim - t for t in tps]

        # Get the full y, x coords of the origin, and window coords of date
        if y in y_edges_1:
            ycenter = y
            oy = 0
        elif y in y_edges_2:
            ycenter = y - ydim
            oy = y - self.spatial_param
        else:
            ycenter = self.spatial_param
            oy = y - self.spatial_param
        if x in x_edges_1:
            xcenter = x
            ox = 0
        elif x in x_edges_2:
            xcenter = x - xdim
            ox = x - self.spatial_param
        else:
            xcenter = self.spatial_param
            ox = x - self.spatial_param
        center = [ycenter, xcenter]
        origin = [oy, ox]
        return top, bottom, left, right, center, origin

    def get_availables(self):
        """
        To save time, avoid checking cells with no events at any time step.
        Create a mask of max values at each point. If the maximum at a cell is
        less than or equal to zero there were no values and it will not be
        checked in the event classification step.
        """
        # Low memory - Somehow leads to slow loop in get_event_perimeters_3d
#        # We want to get the mask without pulling the whole thing into memory
#        burns = xr.open_dataset(self.nc_path, chunks={"x": 500, "y": 500})
#
#        # Pull in only the single max value array
#        mask = burns.max(dim="time").compute()
#
#        # Get the y, x positions where one or more burns were detected
#        locs = np.where(mask.value.values > 0)

#        # Now pair these
#        available_pairs = []
#        for i in range(len(locs[0])):
#            available_pairs.append([locs[0][i], locs[1][i]])

#        # Leaving the data set open causes problems
#        burns.close()

        # Using memory - can handle large tiles, but gets pretty high
        mask = self.input_array.max(dim="time")
        locs = np.where(mask > 0)
        available_pairs = []
        for i in range(len(locs[0])):
            available_pairs.append([locs[0][i], locs[1][i]])
        del mask

        # Done.
        return available_pairs

    def get_event_perimeters_3d(self):
        """
        Iterate through each cell in the 3D MODIS Burn Date tile and group it
        into fire events using the space-time window.

        self = EventGrid("/home/travis/github/firedpy/data/rasters/burn_area/netcdfs/h21v08.nc")
        """
        print("Filtering out cells with no events...")
        available_pairs = self.get_availables()

        # If the other pointers aren"t closed first, this will be very slow
        arr = self.input_array

        # This is to check the window positions
        nz, ny, nx = arr.shape
        dims = [ny, nx]
        perimeters = []

        # traverse spatially, processing each burn day
        print("Building event perimeters...")
        for pair in tqdm(available_pairs, position=0):
            # Separate coordinates
            y, x = pair

            # get the spatial window
            [top, bottom, left,
             right, center, origin] = self.get_spatial_window(y, x, dims)
            cy, cx = center

            # what if we pull in the window?
            window = arr[:, top:bottom+1, left:right+1].data

            # The center of the window is the target burn day
            center_burn = window[:, cy, cx]
            center_burn = center_burn[center_burn > 0]

            # Loop through each event in the window and identify neighbors
            for burn in center_burn:
                new_pts = []
                curr_event_ids = []

                # Now we can get the values and position right away
                diff = abs(burn - window)
                val_locs = np.where(diff <= self.temporal_param)
                y_locs = val_locs[1]
                x_locs = val_locs[2]
                oy, ox = origin

                # Get the actual x, y positions from the window coordinates
                vals = window[val_locs]
                ys = [oy + yl for yl in y_locs]
                xs = [ox + xl for xl in x_locs]

                # Now get the geographic coordinates from tile positions
                all_ys = self.coordinates["y"].data
                all_xs = self.coordinates["x"].data
                ys = all_ys[ys]
                xs = all_xs[xs]

                # Now check if this point is in the event_grid yet
                for i in range(len(vals)):
                    curr_pt = (float(ys[i]), float(xs[i]), float(vals[i]))

                    # already assigned to an event
                    if (curr_pt in self.event_grid):
                        if self.event_grid[curr_pt] not in curr_event_ids:
                            curr_event_ids.append(self.event_grid[curr_pt])
                    else:
                        new_pts.append(curr_pt)

                # If this is a new event
                if len(curr_event_ids)==0:
                    # create a new perimeter object
                    perimeter = EventPerimeter(self.next_event_id, new_pts)

                    # append to perimeters list
                    perimeters.append(perimeter)

                    # add points to the grid
                    self.add_event_grid(self.next_event_id, new_pts)

                    # increment the event ID
                    self.next_event_id += 1

                # If all points part of same existing event
                elif len(curr_event_ids) == 1:
                    event_id = curr_event_ids[0]
                    if len(new_pts):
                        perimeters[event_id - 1].add_coordinates(new_pts)
                        self.add_event_grid(event_id, new_pts)

                # events overlap
                else:
                    perimeters = self.merge_perimeters(perimeters,
                                                       curr_event_ids[0],
                                                       curr_event_ids[1])

        return perimeters
