-------------------

ABSTRACT

-------------------

This is event- and daily-level polygons for the Fire event delineation (FIRED) product for ITALY from November 2000  to July 2021. It is derived from the MODIS MCD64A1 burned area product (see https://lpdaac.usgs.gov/products/mcd64a1v006/ for more details). The MCD64A1 is a monthly raster grid of estimated burned dates. Firedpy (www.github.com/earthlab/firedpy) is an algorithm that converts these rasters into events by stacking the entire time series into a spatial-temporal data cube, then uses an algorithm to assign event identification numbers to pixels that fit into the same 3-dimensional spatial temporal window. This particular dataset was created using a spatial parameter of 5 pixels and 11 days. The primary benefit to this dataset over others is the ability to calculate fire spread rate. For each of these products (events and daily) the event identification numbers are the same, but the event-level product has only single polygons for each entire event, while the daily product has separate polygons for each date per event.  See the accompanying metadata files for the statistics provided by each data set. See the associated paper for more details on the methods and more:

Balch, J.K.; St. Denis, L.A.; Mahood, A.L.; Mietkiewicz, N.P.; Williams, T.M.; McGlinchy, J.; Cook, M.C. FIRED (Fire Events Delineation): An Open, Flexible Algorithm and Database of US Fire Events Derived from the MODIS Burned Area Product (2001–2019). Remote Sens. 2020, 12, 3498. https://doi.org/10.3390/rs12213498 

-------------------

GENERAL INFORMATION

-------------------



1. Title of Dataset:  FIRED-  ITALY 



2. Authors: Jennifer K. Balch, Lise A. St. Denis, Adam L. Mahood, Nathan P.  Mietkiewicz, Travis Williams, Joe McGlinchy, Maxwell C. Cook, Estelle J. Lindrooth.



3. Contact information: jennifer.balch@colorado.edu; adam.mahood@colorado.edu



4. Date of data collection:November 2000 - July 2021



--------------------------

SHARING/ACCESS INFORMATION

--------------------------



1. Licenses/restrictions placed on the data: MIT



2. Links to publications that cite or use the data: TBD



3. Links to other publicly accessible locations of the data: None



4. Recommended citation for the data: 


Balch, J.K.; St. Denis, L.A.; Mahood, A.L.; Mietkiewicz, N.P.; Williams, T.M.; McGlinchy, J.; Cook, M.C. FIRED (Fire Events Delineation): An Open, Flexible Algorithm and Database of US Fire Events Derived from the MODIS Burned Area Product (2001–2019). Remote Sens. 2020, 12, 3498. https://doi.org/10.3390/rs12213498
 
-------------------

DATA & FILE OVERVIEW

-------------------

1. File List: 

     1. Tables: 
         A. fired_italy_to2021182_events.csv

         B. fired_italy_to2021182_daily.csv

             i. This is each fired event split into daily polygons. Each polygon will have an id for the event (which may encompass multiple polygons), and a unique date.

     2. Shapefiles: 
 
         A. fired_italy_to2021182_events.gpkg

         B. fired_italy_to2021182_daily.gpkg

             i. This is each fired event split into daily polygons. Each polygon will have an id for the event (which may encompass multiple polygons), and a unique date.

-------------------

METHODOLOGICAL INFORMATION

-------------------

 1. Spatial window: 1 

 2. Temporal window: 5 

See Balch et al 2020 for complete methods. DOI: https://doi.org/10.3390/rs12213498

-------------------

DATA-SPECIFIC INFORMATION FOR: fired_italy_to2021182_events.csv and fired_italy_to2021182_events.gpkg
-------------------

1. Number of variables: 24

2. Number of cases/rows: 10230

3. Projection information (proj4 string): +proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs

3.1. The projection is the native projection from the MODIS MCD64A1 burned area product from which this dataset is derived. The MCD64A1 product is a raster grid with a resolution of 463 meters. More info at https://lpdaac.usgs.gov/products/mcd64a1v006/ 

4. Variable List: 

        A. Name: id  

        	i. Description: Unique identifier of the fire event. 

        B. Name: ig_date 

            i. Description: The earliest date contained in the event  

        C. Name: ig_day  

            i. Description: The day of the year of the earliest date contained in the event  

        D. Name: ig_month  

            i. Description: The month of the earliest date contained in the event  

        E. Name: ig_year 

            i. Description: The year of the earliest date contained in the event.  

        F. Name: last_date  

            i. Description: The latest date contained in the event  

        G. Name: event_day   

            i. Description: Days since ignition date + 1 (ignition date is day 1)   

        H. Name: pixels   

            i. Description: Total number of pixels burned that day.   

        I. Name: tot_px   

            i. Description:  Total pixels burned for the entire event.   

        J. Name: tot_ar_km2   

            i. Description: Area burned in square kilometers for the entire event.   

        K. Name: fsr_px_dy   

            i. Description: Total pixels burned for the entire event divided by the duration of the fire event.   

        L. Name: fsr_km2_dy   

            i. Description: Total kilometers burned for the entire event divided by the duration of the fire event.   

        M. Name: mx_grw_px   

            i. Description: maximum growth in pixels   

        N. Name: mn_grw_px   

            i. Description: minimum growth in pixels   

        O. Name: mu_grw_px   

            i. Description: mean growth in pixels   

        P. Name: mx_grw_km2   

            i. Description: maximum growth in square kilometers   

        Q. Name: mn_grw_km2   

            i. Description: minimum growth in square kilometers   

        R. Name: mu_grw_km2  

            i. Description: mean growth in square kilometers   

        S. Name: mx_grw_dte   

            i. Description: date of maximum   

        T. Name: lc_code   

            i. Description: Numeric code for the landcover type extracted from the MODIS landcover product for the year preceding the fire.  
 
        U. Name: lc_mode   

            i. Description: Numeric code for the landcover type extracted from the MODIS landcover product for the year preceding the fire.   

        V. Name: lc_name   

            i. Description: Character string of the landcover type from the year before the fire.   

        W. Name: lc_desc   

            i. Description: Character string description of the landcover type from the year before the fire.   

        X. Name: lc_type   

            i. Description: Which landcover classification type was used from the MCD12Q1 product? Default is IGBP global vegetation classification scheme   

        Y. Name: eco_mode  

            i. Description: Modal ecoregion code  

        Z. Name: eco_type   

            i. Description: Which type and level of ecoregion classification was used (North america EPA (levels 1-3) vs World Wildlife Federation)   

        AA. Name: eco_name  

            i. Description: Character string of the ecoregion type where the event occurred.   

        BB. Name: ig_utm_x  

            i. Description: estimated ignition x coordinate   

        CC. Name: ig_utm_y  

            i. Description: estimated ignition y coordinate   

        DD. Name: tot_perim   

            i. Description: Total perimeter of the fire event. 

-------------------

DATA-SPECIFIC INFORMATION FOR: fired_italy_to2021182_daily.csv and fired_italy_to2021182_daily.gpkg
-------------------

1. Number of variables: 29

2. Number of cases/rows: 24202

3. Projection information (proj4 string): +proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs

3.1. The projection is the native projection from the MODIS MCD64A1 burned area product from which this dataset is derived. The MCD64A1 product is a raster grid with a resolution of 463 meters. More info at https://lpdaac.usgs.gov/products/mcd64a1v006/ 

4. Variable list:  

            EE. Name: id 

    		      i. Description: Unique identifier of the fire event. 

    		FF. Name: did  

    		      i. Description: unique identifier of the day within the event  

    		GG. Name: date   

    		      i. Description: date that the area burned  

    		HH. Name: ig_date  

    		      i. Description: The earliest date contained in the event  

    		II. Name: ig_day   

    		      i. Description: The day of the year of the earliest date contained in the event  

    		JJ. Name: ig_month   

    		      i. Description: The month of the earliest date contained in the event  

    		KK. Name: ig_year  

    		      i. Description: The year of the earliest date contained in the event.  

    		LL. Name: last_date  

    	           i. Description: The latest date contained in the event  

    		MM. Name: event_day   

    		      i. Description: Days since ignition date + 1 (ignition date is day 1)  

    		NN. Name: pixels  

    		      i. Description: Total number of pixels burned that day.   

    		OO. Name: tot_px   

    		      i. Description:  Total pixels burned for the entire event.   

    		PP. Name: dy_ar_km2 -  

    		      i. Description:  Area burned in square kilometers that day.   

    		QQ. Name: tot_ar_km2  

    		      i. Description: Area burned in square kilometers for the entire event.  
 
    		RR. Name: fsr_px_dy  

    		      i. Description: Total pixels burned for the entire event divided by the duration of the fire event.   

    		SS. Name: fsr_km2_dy  

    		      i. Description: Total kilometers burned for the entire event divided by the duration of the fire event.   

    		TT. Name: mx_grw_px  

    		      i. Description: Maximum daily fire growth per event in pixels  

    		UU. Name: mn_grw_px  

    		      i. Description: Minimum daily fire growth per event in pixels  

    		VV. Name: mu_grw_px  

    		      i. Description: Mean daily fire growth per event in pixels  

    		WW. Name: mx_grw_km2  

    		      i. Description: Maximum daily fire growth per event in square kilometers  

    		XX. Name: mn_grw_km2  

    		      i. Description: Minimum daily fire growth per event in square kilometers  

    		YY. Name: mu_grw_km2  

    		      i. Description: Mean daily fire growth per event in square kilometers  

    		ZZ. Name: mx_grw_dte   

    		      i. Description: Date of maximum fire growth  

    		AAA. Name: lc_code  

    		      i. Description: Numeric code for the landcover type extracted from the MODIS landcover product for the year preceding the fire.   

    		BBB. Name: lc_mode  

    		      i. Description: Numeric code for the landcover type extracted from the MODIS landcover product for the year preceding the fire.  

    		CCC. Name: lc_name  

    		      i. Description: Character string of the landcover type from the year before the fire.   

    		DDD. Name: lc_desc  

    		      i. Description: Character string description of the landcover type from the year before the fire.   

    		EEE. Name: lc_type  

    		      i. Description: The landcover classification scheme used  

    		FFF. Name: eco_mode  

    		      i. Description: modal ecoregion type  

    		GGG. Name: eco_type  

    		      i. Description: modal ecoregion type  

    		HHH. Name: eco_name  

    		      i. Description: Character string of the landcover type from the year before the fire.   

    		III. Name: ig_utm_x  

    		      i. Description: estimated ignition x coordinate  

    		JJJ. Name: ig_utm_y  

    		      i. Description: estimated ignition y coordinate 
