import csv
import datetime as dt
import re
import os
import pandas as pd

def makeReadMe(proj_dir, tilename, file_base, input, first_date, last_date, ecoregion_type, ecoregion_level, landcover_type, daily, spatial_param, temporal_param, shapefile, shp_type):
    read_path = os.path.join(proj_dir, 'outputs',file_base+"_README.txt")
    last_date = str(last_date)
    first_date = str(first_date)
    file_name = file_base+"_to"+last_date
    last_year = dt.datetime(year=int(last_date[:4]), month=1, day=1)
    last_date = last_year + dt.timedelta(int(last_date[4:]))
    first_year = dt.datetime(year=int(first_date[:4]), month=1, day=1)
    first_date = first_year + dt.timedelta(int(first_date[4:]))
    first_event = first_date.strftime("%B %Y")
    last_event = last_date.strftime("%B %Y")
    if input == 3:
        name = tilename
    else:
        tilename = tilename.upper()
        name = tilename.split("_")

    with open(read_path, "w") as text_file:
        print("-------------------\n", file=text_file)
        print("ABSTRACT\n", file=text_file)
        print("-------------------\n", file=text_file)
        if input == 3:
            print("This is event- and daily-level polygons for the Fire event delineation (FIRED) product for MODIS grid tiles", *name , "from {}  to {}. It is derived from the MODIS MCD64A1 burned area product (see https://lpdaac.usgs.gov/products/mcd64a1v006/ for more details). The MCD64A1 is a monthly raster grid of estimated burned dates. Firedpy (www.github.com/earthlab/firedpy) is an algorithm that converts these rasters into events by stacking the entire time series into a spatial-temporal data cube, then uses an algorithm to assign event identification numbers to pixels that fit into the same 3-dimensional spatial temporal window. This particular dataset was created using a spatial parameter of 5 pixels and 11 days. The primary benefit to this dataset over others is the ability to calculate fire spread rate. For each of these products (events and daily) the event identification numbers are the same, but the event-level product has only single polygons for each entire event, while the daily product has separate polygons for each date per event.  See the accompanying metadata files for the statistics provided by each data set. See the associated paper for more details on the methods and more:\n".format(first_event, last_event), file=text_file)
        else:
            print("This is event- and daily-level polygons for the Fire event delineation (FIRED) product for", *name , "from {}  to {}. It is derived from the MODIS MCD64A1 burned area product (see https://lpdaac.usgs.gov/products/mcd64a1v006/ for more details). The MCD64A1 is a monthly raster grid of estimated burned dates. Firedpy (www.github.com/earthlab/firedpy) is an algorithm that converts these rasters into events by stacking the entire time series into a spatial-temporal data cube, then uses an algorithm to assign event identification numbers to pixels that fit into the same 3-dimensional spatial temporal window. This particular dataset was created using a spatial parameter of 5 pixels and 11 days. The primary benefit to this dataset over others is the ability to calculate fire spread rate. For each of these products (events and daily) the event identification numbers are the same, but the event-level product has only single polygons for each entire event, while the daily product has separate polygons for each date per event.  See the accompanying metadata files for the statistics provided by each data set. See the associated paper for more details on the methods and more:\n".format(first_event, last_event), file=text_file)
        print("Balch, J.K.; St. Denis, L.A.; Mahood, A.L.; Mietkiewicz, N.P.; Williams, T.M.; McGlinchy, J.; Cook, M.C. FIRED (Fire Events Delineation): An Open, Flexible Algorithm and Database of US Fire Events Derived from the MODIS Burned Area Product (2001–2019). Remote Sens. 2020, 12, 3498. https://doi.org/10.3390/rs12213498 \n", file=text_file)
        print("""-------------------\n
GENERAL INFORMATION\n
-------------------\n


1. Title of Dataset:  FIRED- """, *name, """\n


2. Authors: Jennifer K. Balch, Lise A. St. Denis, Adam L. Mahood, Nathan P.  Mietkiewicz, Travis Williams, Joe McGlinchy, Maxwell C. Cook, Estelle J. Lindrooth.\n


3. Contact information: jennifer.balch@colorado.edu; adam.mahood@colorado.edu\n


4. Date of data collection:{} - {}\n


--------------------------\n
SHARING/ACCESS INFORMATION\n
--------------------------\n


1. Licenses/restrictions placed on the data: MIT\n


2. Links to publications that cite or use the data: TBD\n


3. Links to other publicly accessible locations of the data: None\n


4. Recommended citation for the data: \n

Balch, J.K.; St. Denis, L.A.; Mahood, A.L.; Mietkiewicz, N.P.; Williams, T.M.; McGlinchy, J.; Cook, M.C. FIRED (Fire Events Delineation): An Open, Flexible Algorithm and Database of US Fire Events Derived from the MODIS Burned Area Product (2001–2019). Remote Sens. 2020, 12, 3498. https://doi.org/10.3390/rs12213498\n """.format(first_event, last_event), file = text_file)



        print("-------------------\n", file=text_file)
        print("DATA & FILE OVERVIEW\n", file=text_file)
        print("-------------------\n", file=text_file)
        print("1. File List: \n", file=text_file)
        if daily == 'yes' and shapefile:
            print("     1. Tables: ", file=text_file)
            print("         A. {}_events.csv\n".format(file_name), file=text_file)
            print("         B. {}_daily.csv\n".format(file_name), file=text_file)
            print("             i. This is each fired event split into daily polygons. Each polygon will have an id for the event (which may encompass multiple polygons), and a unique date.\n", file=text_file)
            print("     2. Shapefiles: \n ", file=text_file)
            if shp_type == 'gpkg':
                print("         A. {}_events.gpkg\n".format(file_name, last_date), file=text_file)
                print("         B. {}_daily.gpkg\n".format(file_name), file=text_file)
                print("             i. This is each fired event split into daily polygons. Each polygon will have an id for the event (which may encompass multiple polygons), and a unique date.\n", file=text_file)
            if shp_type == 'shp':
                print("         A. {}_events.shp\n".format(file_name), file=text_file)
                print("         B. {}_daily.shp\n".format(file_name), file=text_file)
                print("             i. This is each fired event split into daily polygons. Each polygon will have an id for the event (which may encompass multiple polygons), and a unique date.\n", file=text_file)
            if shp_type == 'both':
                print("         A. {}_events.shp\n".format(file_name), file=text_file)
                print("         B. {}_daily.shp\n".format(file_name), file=text_file)
                print("             i. This is each fired event split into daily polygons. Each polygon will have an id for the event (which may encompass multiple polygons), and a unique date.\n", file=text_file)
                print("         C. {}_events.gpkg\n".format(file_name), file=text_file)
                print("         D. {}_daily.gpkg\n".format(file_name), file=text_file)
                print("             i. This is each fired event split into daily polygons. Each polygon will have an id for the event (which may encompass multiple polygons), and a unique date.\n", file=text_file)

        elif daily == 'yes' and shapefile == False:
            print("     1. Tables:\n ", file=text_file)
            print("         A. {}_events.csv\n".format(file_name), file=text_file)
            print("         B. {}_daily.csv\n".format(file_name), file=text_file)
            print("             i. This is each fired event split into daily polygons. Each polygon will have an id for the event (which may encompass multiple polygons), and a unique date.\n", file=text_file)
        elif daily == 'no' and shapefile:
            print("     1. Table:\n ", file=text_file)
            print("         A. {}_events.csv\n".format(file_name), file=text_file)
            print("     2. Shapefile: \n", file=text_file)
            if shp_type == 'gpkg':
                print("         A. {}_events.gpkg\n".format(file_name), file=text_file)
            if shp_type == 'shp':
                print("         A. {}_events.shp\n".format(file_name), file=text_file)
            if shp_type == 'both':
                print("         A. {}_events.shp\n".format(file_name), file=text_file)
                print("         B. {}_events.gpkg\n".format(file_name), file=text_file)

        else:
            print("     1. Table:\n ", file=text_file)
            print("         A. {}_events.csv\n".format(file_name), file=text_file)
        print("-------------------\n", file=text_file)
        print("METHODOLOGICAL INFORMATION\n", file=text_file)
        print("-------------------\n", file=text_file)
        print(" 1. Spatial window: {} \n".format(spatial_param), file=text_file)
        print(" 2. Temporal window: {} \n".format(temporal_param), file=text_file)
        print("See Balch et al 2020 for complete methods. DOI: https://doi.org/10.3390/rs12213498\n", file=text_file)
        print("-------------------\n", file=text_file)
        if shapefile:
            if shp_type == 'gpkg':
                print("DATA-SPECIFIC INFORMATION FOR: {}_events.csv and {}_events.gpkg".format(file_name, file_name), file=text_file)
            if shp_type == 'shp':
                print("DATA-SPECIFIC INFORMATION FOR: {}_events.csv and {}_events.shp".format(file_name, file_name), file=text_file)
            if shp_type == 'both':
                print("DATA-SPECIFIC INFORMATION FOR: {}_events.csv, {}_events.gpkg, {}_events.shp".format(file_name, file_name,file_name), file=text_file)
        else:
            print("DATA-SPECIFIC INFORMATION FOR: {}_events.csv".format(file_name), file=text_file)

        print("-------------------\n", file=text_file)
        filepath = os.path.join(proj_dir, 'outputs', 'tables',
                                        file_name+"_events.csv")
        reader = pd.read_csv(filepath)
        row_count = len(reader)
        print("1. Number of variables: 24\n", file=text_file)
        print("2. Number of cases/rows: {}\n".format(row_count), file=text_file)
        print("3. Projection information (proj4 string): +proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs\n", file = text_file)
        print("3.1. The projection is the native projection from the MODIS MCD64A1 burned area product from which this dataset is derived. The MCD64A1 product is a raster grid with a resolution of 463 meters. More info at https://lpdaac.usgs.gov/products/mcd64a1v006/ \n", file = text_file)
        print("""4. Variable List: \n
        A. Name: id  \n
        	i. Description: Unique identifier of the fire event. \n
        B. Name: ig_date \n
            i. Description: The earliest date contained in the event  \n
        C. Name: ig_day  \n
            i. Description: The day of the year of the earliest date contained in the event  \n
        D. Name: ig_month  \n
            i. Description: The month of the earliest date contained in the event  \n
        E. Name: ig_year \n
            i. Description: The year of the earliest date contained in the event.  \n
        F. Name: last_date  \n
            i. Description: The latest date contained in the event  \n
        G. Name: event_day   \n
            i. Description: Days since ignition date + 1 (ignition date is day 1)   \n
        H. Name: pixels   \n
            i. Description: Total number of pixels burned that day.   \n
        I. Name: tot_px   \n
            i. Description:  Total pixels burned for the entire event.   \n
        J. Name: tot_ar_km2   \n
            i. Description: Area burned in square kilometers for the entire event.   \n
        K. Name: fsr_px_dy   \n
            i. Description: Total pixels burned for the entire event divided by the duration of the fire event.   \n
        L. Name: fsr_km2_dy   \n
            i. Description: Total kilometers burned for the entire event divided by the duration of the fire event.   \n
        M. Name: mx_grw_px   \n
            i. Description: maximum growth in pixels   \n
        N. Name: mn_grw_px   \n
            i. Description: minimum growth in pixels   \n
        O. Name: mu_grw_px   \n
            i. Description: mean growth in pixels   \n
        P. Name: mx_grw_km2   \n
            i. Description: maximum growth in square kilometers   \n
        Q. Name: mn_grw_km2   \n
            i. Description: minimum growth in square kilometers   \n
        R. Name: mu_grw_km2  \n
            i. Description: mean growth in square kilometers   \n
        S. Name: mx_grw_dte   \n
            i. Description: date of maximum   \n
        T. Name: lc_code   \n
            i. Description: Numeric code for the landcover type extracted from the MODIS landcover product for the year preceding the fire.  \n 
        U. Name: lc_mode   \n
            i. Description: Numeric code for the landcover type extracted from the MODIS landcover product for the year preceding the fire.   \n
        V. Name: lc_name   \n
            i. Description: Character string of the landcover type from the year before the fire.   \n
        W. Name: lc_desc   \n
            i. Description: Character string description of the landcover type from the year before the fire.   \n
        X. Name: lc_type   \n
            i. Description: Which landcover classification type was used from the MCD12Q1 product? Default is IGBP global vegetation classification scheme   \n
        Y. Name: eco_mode  \n
            i. Description: Modal ecoregion code  \n
        Z. Name: eco_type   \n
            i. Description: Which type and level of ecoregion classification was used (North america EPA (levels 1-3) vs World Wildlife Federation)   \n
        AA. Name: eco_name  \n
            i. Description: Character string of the ecoregion type where the event occurred.   \n
        BB. Name: ig_utm_x  \n
            i. Description: estimated ignition x coordinate   \n
        CC. Name: ig_utm_y  \n
            i. Description: estimated ignition y coordinate   \n
        DD. Name: tot_perim   \n
            i. Description: Total perimeter of the fire event. \n""", file = text_file)
        if daily == 'yes':
            print("-------------------\n", file=text_file)
            if shapefile:
                if shp_type == 'gpkg':
                    print("DATA-SPECIFIC INFORMATION FOR: {}_daily.csv and {}_daily.gpkg".format(file_name, file_name), file=text_file)
                if shp_type == 'shp':
                    print("DATA-SPECIFIC INFORMATION FOR: {}_daily.csv and {}_daily.shp".format(file_name, file_name), file=text_file)
                if shp_type == 'both':
                    print("DATA-SPECIFIC INFORMATION FOR: {}_daily.csv, {}_daily.gpkg, {}_daily.shp".format(file_name, file_name,file_name), file=text_file)
            else:
                print("DATA-SPECIFIC INFORMATION FOR: {}_daily.csv".format(file_name), file=text_file)
            print("-------------------\n", file=text_file)
            filepath = os.path.join(proj_dir, 'outputs', 'tables',
                                            file_name+"_daily.csv")
            reader = pd.read_csv(filepath)
            row_count = len(reader)
            print("1. Number of variables: 29\n", file=text_file)
            print("2. Number of cases/rows: {}\n".format(row_count), file=text_file)
            print("3. Projection information (proj4 string): +proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs\n", file = text_file)
            print("3.1. The projection is the native projection from the MODIS MCD64A1 burned area product from which this dataset is derived. The MCD64A1 product is a raster grid with a resolution of 463 meters. More info at https://lpdaac.usgs.gov/products/mcd64a1v006/ \n", file = text_file)
            print("""4. Variable list:  \n
            EE. Name: id \n
    		      i. Description: Unique identifier of the fire event. \n
    		FF. Name: did  \n
    		      i. Description: unique identifier of the day within the event  \n
    		GG. Name: date   \n
    		      i. Description: date that the area burned  \n
    		HH. Name: ig_date  \n
    		      i. Description: The earliest date contained in the event  \n
    		II. Name: ig_day   \n
    		      i. Description: The day of the year of the earliest date contained in the event  \n
    		JJ. Name: ig_month   \n
    		      i. Description: The month of the earliest date contained in the event  \n
    		KK. Name: ig_year  \n
    		      i. Description: The year of the earliest date contained in the event.  \n
    		LL. Name: last_date  \n
    	           i. Description: The latest date contained in the event  \n
    		MM. Name: event_day   \n
    		      i. Description: Days since ignition date + 1 (ignition date is day 1)  \n
    		NN. Name: pixels  \n
    		      i. Description: Total number of pixels burned that day.   \n
    		OO. Name: tot_px   \n
    		      i. Description:  Total pixels burned for the entire event.   \n
    		PP. Name: dy_ar_km2 -  \n
    		      i. Description:  Area burned in square kilometers that day.   \n
    		QQ. Name: tot_ar_km2  \n
    		      i. Description: Area burned in square kilometers for the entire event.  \n 
    		RR. Name: fsr_px_dy  \n
    		      i. Description: Total pixels burned for the entire event divided by the duration of the fire event.   \n
    		SS. Name: fsr_km2_dy  \n
    		      i. Description: Total kilometers burned for the entire event divided by the duration of the fire event.   \n
    		TT. Name: mx_grw_px  \n
    		      i. Description: Maximum daily fire growth per event in pixels  \n
    		UU. Name: mn_grw_px  \n
    		      i. Description: Minimum daily fire growth per event in pixels  \n
    		VV. Name: mu_grw_px  \n
    		      i. Description: Mean daily fire growth per event in pixels  \n
    		WW. Name: mx_grw_km2  \n
    		      i. Description: Maximum daily fire growth per event in square kilometers  \n
    		XX. Name: mn_grw_km2  \n
    		      i. Description: Minimum daily fire growth per event in square kilometers  \n
    		YY. Name: mu_grw_km2  \n
    		      i. Description: Mean daily fire growth per event in square kilometers  \n
    		ZZ. Name: mx_grw_dte   \n
    		      i. Description: Date of maximum fire growth  \n
    		AAA. Name: lc_code  \n
    		      i. Description: Numeric code for the landcover type extracted from the MODIS landcover product for the year preceding the fire.   \n
    		BBB. Name: lc_mode  \n
    		      i. Description: Numeric code for the landcover type extracted from the MODIS landcover product for the year preceding the fire.  \n
    		CCC. Name: lc_name  \n
    		      i. Description: Character string of the landcover type from the year before the fire.   \n
    		DDD. Name: lc_desc  \n
    		      i. Description: Character string description of the landcover type from the year before the fire.   \n
    		EEE. Name: lc_type  \n
    		      i. Description: The landcover classification scheme used  \n
    		FFF. Name: eco_mode  \n
    		      i. Description: modal ecoregion type  \n
    		GGG. Name: eco_type  \n
    		      i. Description: modal ecoregion type  \n
    		HHH. Name: eco_name  \n
    		      i. Description: Character string of the landcover type from the year before the fire.   \n
    		III. Name: ig_utm_x  \n
    		      i. Description: estimated ignition x coordinate  \n
    		JJJ. Name: ig_utm_y  \n
    		      i. Description: estimated ignition y coordinate """, file = text_file)
