# -*- coding: utf-8 -*-
"""
Translate the data frame of wildfire events into polygons.

Notes:
    - Some of the vector manipulations in this step are quite slow. Keep an
      eye out for a stable release of dask-geopandas sometime in late 2019:

          https://github.com/mrocklin/dask-geopandas

Created on Thu Jun 20 09:40:59 2019
@author: Travis
"""
import argparse
import geopandas as gpd
from netCDF4 import Dataset
import pandas as pd
from shapely.geometry import Point, Polygon, MultiPolygon
import shapely.speedups
import time
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

pd.options.mode.chained_assignment = None
shapely.speedups.enable()

parser = argparse.ArgumentParser()
parser.add_argument('-dest1', dest='dest1',
                    default='data/shapefiles/daily_polygons.gpkg',
                    help=("""Use this to specify the location and filename of
                          the daily level shapefile. Defaults to 
                          data/shapefiles/daily_polygons.gpkg"""))
parser.add_argument('-dest2', dest='dest2',
                    default='data/shapefiles/event_polygons.gpkg',
                    help=("""Use this to specify the location and filename of
                          the event level shapefile. Defaults to
                          data/shapefiles/event_polygons.gpkg"""))
parser.add_argument('-src', dest='src',
                    default='data/tables/modis_events.csv',
                    help=("""Specify the location and filename of the fire
                          events data frame. Defaults to
                          data/tables/modis_events.csv"""))
args = parser.parse_args()
dest1 = args.dest1
dest2 = args.dest2
src = args.src

# Start the timer (seconds)
start = time.perf_counter()

# Read in the event table and a reference nc file for geometric information
print('Reading classified fire event table...')
df = pd.read_csv(src)

# Go ahead and create daily id (did) for later
df['did'] = df['id'].astype(str) + '-' + df['date']

# Use a sample for geometries
sample = Dataset('data/rasters/burn_area/netcdfs/h08v04.nc')

#conus = gpd.read_file('data/shapefiles/conus.shp')
#conus = conus.rename(columns={'NAME': 'state'})
crs = sample.variables['crs']
geom = crs.geo_transform
proj4 = crs.proj4
res = [geom[1], geom[-1]]

# Filter columns, center pixel coordinates, and remove repeating pixels
df = df[['id', 'did', 'date', 'x', 'y']]
df['x'] = df['x'] + (res[0]/2)
df['y'] = df['y'] + (res[1]/2)

# Each entry in the df gets a point object from the x and y coordinates.
print('Converting data frame to spatial object...')
df['geometry'] = df[['x', 'y']].apply(lambda x: Point(tuple(x)), axis=1)
gdf = df[['id', 'did', 'date', 'geometry']]
gdf = gpd.GeoDataFrame(gdf, crs=proj4, geometry=gdf['geometry'])

# Create a circle buffer
print("Creating buffer...")
geometry = gdf.buffer(1 + (res[0]/2))
gdf['geometry'] = geometry

# Then create a square envelope around the circle
gdf['geometry'] = gdf.envelope

# Now add the first date of each event and merge daily event detections
print("Dissolving polygons...")
gdf['start_date'] = gdf.groupby('id')['date'].transform('min')  
gdfd = gdf.dissolve(by='did', as_index=False)
gdfd['year'] = gdfd['start_date'].apply(lambda x: x[:4])
gdfd['month'] = gdfd['start_date'].apply(lambda x: x[5:7])
#gdfd = gdfd[gdfd['month'].isin(['06', '07', '08', '09'])]

# Save the daily before dissolving into event level
print('Saving daily file...')
gdfd.to_file(dest1, driver='GPKG')

# Now merge into event level polygons
gdf = gdf[['id', 'start_date', 'geometry']]
gdf = gdf.dissolve(by='id', as_index=False)

# For each geometry, if it is a single polygon, cast as a multipolygon
print("Converting polygons to multipolygons...")
def asMultiPolygon(polygon):
    if type(polygon) == Polygon:
        polygon = MultiPolygon([polygon])
    return polygon
gdf['geometry'] = gdf['geometry'].apply(asMultiPolygon)

# Calculate perimeter length
print('Calculating perimeter lengths...')
gdf['final_perimeter'] = gdf['geometry'].length  # <--------------------------- Check accuracy of this, QGIS is slightly different (also check adams)

# Now save as a geopackage  # <------------------------------------------------ Should we also make a shapefile for ESRI users? Make it optional.
print('Saving event-level file...')
gdf.to_file(dest2, driver='GPKG')

# Print the time it took
end = time.perf_counter()
seconds = end - start
minutes = seconds/60
print('Job completed in {} minutes'.format(round(minutes, 2)))
