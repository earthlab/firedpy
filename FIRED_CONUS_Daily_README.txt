-------------------
GENERAL INFORMATION
-------------------


1. Title of Dataset:  FIRED CONUS: Daily


2. Authors: Jennifer K. Balch, Lise A. St. Denis, Adam L. Mahood, Nathan P.  Mietkiewicz, Travis Williams, Joe McGlinchy, Maxwell C. Cook.


3. Contact information: jennifer.balch@colorado.edu; adam.mahood@colorado.edu


4. Date of data collection: 2001 - January 2019


--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 


1. Licenses/restrictions placed on the data: MIT


2. Links to publications that cite or use the data: TBD


3. Links to other publicly accessible locations of the data: None


4. Recommended citation for the data: TBD


---------------------
DATA & FILE OVERVIEW
---------------------


1. File List:


   A. Filename: fired_conus_daily_nov2001-jan2019.gpkg
	i. This is each fired event split into daily polygons. Each polygon will have an id for the event (which may encompass multiple polygons), and a unique date. 



--------------------------
METHODOLOGICAL INFORMATION
--------------------------

See Balch et al 2020 preprint. DOI: https://doi.org/10.32942/osf.io/nkzpg


-----------------------------------------
DATA-SPECIFIC INFORMATION FOR: fired_conus_daily_nov2001-jan2019.gpkg
-----------------------------------------

1. Number of variables: 24

2. Number of cases/rows: 281,615

3. Projection information (proj4 string): "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

3.1. The projection is the native projetion from the MODIS MCD64A1 burned area product from which this dataset is derived. The MCD64A1 product is a raster grid with a resultion of 463 meters. More info at https://lpdaac.usgs.gov/products/mcd64a1v006/

4. Variable List:
	A. Name: id
	   Description: Unique identifier of the fire event 

        B. Name: date
       	   Description: The date that the polygon burned

	C. Name: pixels
	   Description: Total number of pixels burned that day

	D. Name: l1_eco
	   Description: numeric code for the level 1 ecoregion

	E. Name: lc 
	   Description: numeric code for the landcover type extracted from the MODIS landcover product for the year preceeding the fire

	F. Name: cum_pixels
	   Description: cumulative pixels for the event up to and including that date

	G. Name: total_pixels
	   Description: total pixels burned for the entire event

        H. Name: ignition_date
       	   Description: The earliest date contained in the event

	I. Name: last_date
	   Description: The latest date contained in the event

	J. Name: duration
	   Description: last_date - ignition_date + 1; duration of the fire event

	K. Name: simple_fsr_pixels
	   Description: Total pixels burned for the entire event divided by the duration of the fire event

        L. Name: simple_fsr_km2
       	   Description: Total area burned for the entire event in square kilometers divided by the duration of the fire event

	M. Name: daily_area_km2
	   Description: Area burned in square kilometers that day

        N. Name: cum_area_km2
       	   Description: Cumulative area burned in kilometers squared up to and including that day

	O. Name: total_area_km2
	   Description: Total area burned for the entire event in square kilometers

        P. Name: pct_total_area
       	   Description: Area burned that day divided by the total area of the event, times 100

	Q. Name: pct_cum_area
	   Description: Area burned that day divided by the cumulative area up to including that day, times 100

        R. Name: event day
       	   Description: Days since ignition date + 1 (ignition date is day 1)

	S. Name: ratio_area_added_to_average
	   Description: Area burned that day divided by simple fsr

        T. Name: prior_pixels
       	   Description: How many pixels burned prior to this date

	U. Name: rel_fsr_per_day
	   Description: How many pixels burned on the day divided by the prior pixels
        
	V. Name: l1_ecoregion
       	   Description: Character string of the level 1 ecoregion
        
	W. Name: lc_name
       	   Description: Character string of the landcover type from the year before the fire

	X. Name: geom
	   Description: The geometry information