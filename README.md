

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5507688.svg)](https://doi.org/10.5281/zenodo.5507688)



[![Docker Automated build](https://img.shields.io/docker/automated/earthlab/firedpy?style=plastic)](https://hub.docker.com/repository/docker/earthlab/firedpy/builds)
![GitHub contributors](https://img.shields.io/github/contributors/earthlab/firedpy)
[![GitHub issues](https://img.shields.io/github/issues/earthlab/firedpy)](https://github.com/earthlab/firedpy/issues)
![GitHub commit activity](https://img.shields.io/github/commit-activity/w/earthlab/firedpy) 
![pytests](https://github.com/earthlab/firedpy/workflows/Pytests/badge.svg)
[![Coverage](https://codecov.io/github/earthlab/firedpy/coverage.svg?branch=main)](https://codecov.io/gh/earthlab/firedpy)

# FIREDpy - FIRe Event Delineation for Python

A Python Command-Line Interface for classifying fire events from the Collection 6 MODIS Burned Area Product.

This package uses a space-time window to classify individual burn detections from late 2001 to near-present into discrete events and return both a data table and shapefiles of these events. The user is able to specify the spatial and temporal parameters of the window, as well as the spatial and temporal extent, using either a shapefile or a list of MODIS Sinusoidal Projection tile IDs. Countries, continents and US states are included. Any area from the world may be selected. However, in the current version, memory constraints may limit the extent available for a single model run. Equatorial regions have much more fire activity, and may require much more RAM to process than a normal laptop will have.

The algorithm outputs shapefiles of delineated fire events in either .shp or .gpkg format. In addition to the full event polygons created by default, and the user has the option of having firedpy produce daily-level perimeters, providing a representation of both final and expanding event perimeters. 

<img width="1900" height="544" alt="image" src="https://github.com/user-attachments/assets/0144541b-6fc7-4718-b620-4347ae77881f" />

*Illustration of the event-level and daily-level output of FIREDpy for the 2013 Rim Fire in California. Figure is from Mahood et al. 2022.*


## The MODIS Sinusioidal Grid

In general, the area of intrest for which FIREDpy creates fire perimter products is specified with a shapefile. However, one may also specify study areas for FIREDpy using individual MODIS tiles. This product uses a Sinusoidal, Lambert Azimuthal Equal-Area projection. For more information about this projection, see https://modis-land.gsfc.nasa.gov/MODLAND_grid.html. If you would like to use the MODIS sinusoidal grid IDs to identify FIREDpy study areas, see the visualization of this grid below.

<img width="1900" alt="image" src="https://github.com/earthlab/firedpy/blob/tw/installation/firedpy/data/images/modis_land_id_map_robinson.png" />



## FIREDpy citations

### Methodological information for FIREDpy 1.0:

Balch, J. K., St. Denis, L. A., Mahood, A. L., Mietkiewicz, N. P., Williams, T. W., McGlinchy J, and Cook, M. C. 2020. FIRED (Fire Events Delineation): An open, flexible algorithm & database of U.S. fire events derived from the MODIS burned area product (2001-19). Remote Sensing, 12(21), 3498; https://doi.org/10.3390/rs12213498

### Description of the 2000 - 2021 data sets: 

Mahood, A.L. Lindrooth, E.J., Cook, M.C. and Balch, J.K. 2022. Country-level fire perimeter datasets (2001-2021). Nature Scientific Data, 9(458). https://doi.org/10.1038/s41597-022-01572-3

### Methodological information for FIREDpy 2.0, description of 2000-2025 datasets:

Coming soon...

## Data Sharing Agreement
FIREDpy is currently in active development, and newer versions of the algorithm and data products are shared on an individual basis. All data products generated from unpublished versions of FIREDpy require permission from PI Jennifer K. Balch prior to use in publications, presentations, or public dissemination. These conversations ensure appropriate acknowledgment of the development team's contributions and proper context for the algorithm's current capabilities and limitations. Please use the above citation for attributing credit. For data requests or collaboration inquiries, please contact Nate Hofford (nate.hofford@colorado.edu), University of Colorado Boulder.

## Changes
 - 10/14/2024 FIREDpy V2.0
    - No longer using setup.py. See new instructions below for running it with Docker or installing it locally.
    - Improved fire grouping
    - Improved CLI
    - Access to MODIS burn area product Version 6.1

### BUG ALERT: 

Many of the data products created in Fall 2021 may be shifted by a half pixel, and may lack a coordinate reference system. 

The problem is now fixed, so this will not affect new iterations of firedpy.

Sometimes the server that houses the MCD64A1 product used by firedpy is down. If this happens, you just need to wait until it comes back up.

See the issues tab for more bugs, or to report a new bug!

## Have you used firedpy?

The algorithm and derived data products are under active development. Please raise an issue with any suggestions to help us improve firedpy. Or just email admahood@duck.com and Adam will be overjoyed to talk about firedpy.

## Current status of created products

Already-created products are linked below. They are housed in the CU Scholar data repository in the [Earth Lab Data collection](https://scholar.colorado.edu/collections/pz50gx05h), or [here](https://scholar.colorado.edu/catalog?f%5Bacademic_affiliation_sim%5D%5B%5D=Earth+Lab&locale=en). 

All of the created products have an event-level shapefile in .gpkg and .shp formats. Most countries also have the daily-level shapefile for the V2024 products, but for the V2022 products these were not created for most countries in Africa and Asia due to computational restrictions. 

## Click on a link below to download a fire perimeter dataset for your favorite country

|| **Firedpy version**| 1.0|2.0|2.0|
|-------|---------|----------------------------------------|----------------------------------------|---------------------------|
| **Region** | **Country** | **Version 2022v (variable space-time parameters)** | **Version 2024f (1 pixel 5 days for everything)** | **Version 2025v (variable ST parameters)** |
| **North America** |||||
|| Belize || **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/q524jq679)** ||
|| Canada| [2001-2021](https://scholar.colorado.edu/concern/datasets/gf06g388c) | **[2001-2024](https://scholar.colorado.edu/concern/datasets/0p096866d)** ||
|| Costa Rica || **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/mc87pr85g)** ||
|| El Salvador || **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/m613n034n)**||
|| Guatemala || **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/nv935460v)**||
|| Honduras || **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/n296x082p)**||
|| Mexico || **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5x21th07j)**||
|| Nicaragua || **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/rb68xd57f)**||
|| Panama || **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/6q182m92n)**||
|| USA plus Canada | [2001-2021](https://scholar.colorado.edu/concern/datasets/8336h304x)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/0g354g88p)**||
|| USA: Coterminous + Alaska | [2001-2021](https://scholar.colorado.edu/concern/datasets/d504rm74m) | **[November 2001- July 2024](https://scholar.colorado.edu/concern/datasets/fx719p11c)**||
|| USA: Hawaii | [2001-2021](https://scholar.colorado.edu/concern/datasets/7h149r06p)|**[November 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/sn00b078q)**||
|Mexico and Central America |Belize, Guatemala, Honduras, El Salvador, Nicaragua, Costa Rica, Panama| [2001-2021](https://scholar.colorado.edu/concern/datasets/vd66w1102)||
|**Carribean Sea** |||||
|Whole region |Barbados, Bahamas, Cayman Islands, Cuba, Dominican Republic, Haiti, Jamaica, Montserrat, Puerto Rico, Saint Kitts And Nevis, Trinidad And Tobago, British Virgin Islands, Guadeloupe, Saint Barthelemy| [2001-2021](https://scholar.colorado.edu/concern/datasets/x633f230f) ||
|| Aruba||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/9880vs496)**|| 
||  Antigua and Barbuda||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/st74cs33x)**|| 
 ||  The Bahamas||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/t722hb658)**|| 
 ||  Barbados||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/m613n033c)**|| 
 ||  British Virgin Islands||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/jq085m540)**|| 
 ||  Cayman Islands||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/0r967557n)**|| 
 ||  Cuba||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5712m8365)**|| 
 ||  Curacao||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/ws859h28b)**|| 
 ||  Dominica||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/hd76s170r)**|| 
 ||  Dominican Republic||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1j92g9289)**|| 
 ||  Grenada||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/9880vs50z)**|| 
 ||  Haiti||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/707959237)**|| 
 ||  Jamaica||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/xg94hr218)**|| 
 ||  Montserrat||  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/3r074w69k)**|| 
 ||  USA: Puerto Rico||  **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/gt54kq07d)**|| 
 ||  Saint Barthelemy||  **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/h989r4943)**|| 
 ||  Saint Kitts and Nevis||  **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2f75r969z)**|| 
 ||  Saint Martin||  **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/7d278v61h)**|| 
 ||  Trinidad and Tobago||  **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/9c67wp39g)**|| 
 ||  Turks and Caicos||  **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/cv43nz30q)**|| 
 ||  US Virgin Islands||  **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/pv63g214m)**|| 
|**South America** |||||
||  Argentina| [2001-2021](https://scholar.colorado.edu/concern/datasets/5t34sk58k)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/td96k4315)**||
||  Brazil| [2001-2021](https://scholar.colorado.edu/concern/datasets/05741s90q) |**[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/b8515q25m)**||
||  Bolivia| [2001-2021](https://scholar.colorado.edu/concern/datasets/b2773w83t)|**[November 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/2f75r992r)**||
||  Chile| [2001-2021](https://scholar.colorado.edu/concern/datasets/qr46r2052) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/g732db604)**|| 
||  Colombia| [2001-2021](https://scholar.colorado.edu/concern/datasets/mp48sd91d) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/8c97ks03f)**|| 
|| Ecuador| [2001-2021](https://scholar.colorado.edu/concern/datasets/pc289k34n)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5138jg700)**|| 
||  Guyana|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5h73pz00k)**||
||  Paraguay| [2001-2021](https://scholar.colorado.edu/concern/datasets/rb68xd05p)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/fn1070729)**|| 
||  Peru| [2001-2021](https://scholar.colorado.edu/concern/datasets/x346d5616)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/wh246t78k)**|| 
||  Suriname|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/n583xw60q)**|| 
||  Uruguay| [2001-2021](https://scholar.colorado.edu/concern/datasets/q524jq130) |**[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/q237ht604)**|| 
||  Venezuela| [2001-2021](https://scholar.colorado.edu/concern/datasets/7m01bm95m) |**[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/c534fq76h)**|| 
|  Northern South America |Suriname, French Guiana, Guyana| [2001-2021](https://scholar.colorado.edu/concern/datasets/qv33rx839)||| 
|Entire Western hemisphere, intended for use in conjunction with GOES16 active fire detections. || [Jan 2017 to March 2020](https://scholar.colorado.edu/concern/datasets/d217qq78g)|||
| **Europe** |||||
|| Aland|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5138jg697)**||
|| Albania|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/wp988m53f)**||
|| Andorra|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/xp68kh87x)**||
|| Belarus|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/r781wh632)**||
|| Belgium||**[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/hh63sx489)**||
|| Bosnia and Herzegovina|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/dz010s09q)**||
|| Croatia|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/kh04dr144)**||
|| Cyprus|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/tx31qk48f)**||
|| Czechia|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/fq977w613)**||
|| Estonia|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/73666649m)**||
|| Finland| [2001-2021](https://scholar.colorado.edu/concern/datasets/6395w836j)| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/pk02cc33v)**||
|| Germany|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/6w924d770)**||
|| Greece| [2001-2021](https://scholar.colorado.edu/concern/datasets/bc386k355)| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/z316q349s)**||
|| Greenland|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/bv73c218z)**||
|| Guernsey|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/bz60cx90j)**||
|| Hungary||**[Nov 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/bc386m008)**||
|| Ireland||**[Nov 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/4m90dx459)**||
|| Italy|[2001-2021](https://scholar.colorado.edu/concern/datasets/v979v416g)| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/tb09j7441)**||
|| Kosovo|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/f1881n864)**||
|| Latvia||**[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/fn1070711)**||
|| Liechtenstein|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/f7623f56w)**||
|| Lithuania|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/qb98mg925)**||
|| Macedonia|| **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/g158bj86s)**||
|| Malta|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/sx61dp12d)**||
|| Moldova|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/d217qr276)**||
|| Monaco|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2514nn122)**||
|| Montenegro|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1g05fd436)**||
|| Netherlands|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/0k225c86c)**||
|| Northern Cyprus|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/4m90dx238)**||
|| Norway|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/70795924h)**||
|| Poland| [2001-2021](https://scholar.colorado.edu/concern/datasets/6969z264g)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/4b29b7569)**||
|| Portugal|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/f4752j497)**||
|| Romania|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/w3763834n)**||
|| San Marino|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/zk51vj568)**||
|| Serbia|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5h73pz01v)**||
|| Slovakia|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/cr56n2793)**||
|| Slovenia|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/tq57ns75v)**||
|| Spain|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/xw42n953k)**||
|| Sweden|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1j92g930b)**||
|| Switzerland|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1g05fd42x)**||
|| Ukraine|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/08612q18f)**||
|| United Kingdom|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/sf268684c)**||
|The British Isles|UK and Ireland| [2001-2021](https://scholar.colorado.edu/concern/datasets/pc289k33c)|||
|The Iberian Peninsula| Spain & Portugal| [2001-2021](https://scholar.colorado.edu/concern/datasets/gb19f7006)|||
| Western Europe |France, Germany, Poland, Switzerland, Belgium, Netherlands, Luxembourg and Austria| [2001-2021](https://scholar.colorado.edu/concern/datasets/v692t736f)|||
 | Central to Southern Europe |Estonia, Latvia, Lithuania, Belarus, Ukraine, Czech Republic, Slovakia, Hungary, Romania, Bulgaria, Montenegro, Bosnia, Turkey, Republic Of Moldova, Serbia, Albania, Slovenia, and North Macedonia| [2001-2021](https://scholar.colorado.edu/concern/datasets/7h149r07z)|||
 | Northern Europe |Iceland, Sweden, Norway, and Denmark| [2001-2021](https://scholar.colorado.edu/concern/datasets/sb397945f)|||
| **Africa**|||||
|| Algeria|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/bc386k800)**||
|| Angola| [2001-2021](https://scholar.colorado.edu/concern/datasets/t435gf21z) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/d504rn23b)**||
|| Benin| [2001-2021](https://scholar.colorado.edu/concern/datasets/z603qz58m) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/fb494b180)**||
|| Botswana| [2001-2021](https://scholar.colorado.edu/concern/datasets/b8515p69g)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/nv9354593)**||
|| Burundi| [2001-2021](https://scholar.colorado.edu/concern/datasets/3f462659h) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/vt150k986)**||
|| Burkina Faso| [2001-2021](https://scholar.colorado.edu/concern/datasets/9g54xj875)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/x633f2869)**||
|| Cabo Verde||**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/v692t782j)**||
|| Cameroon| [2001-2021](https://scholar.colorado.edu/concern/datasets/x920fz208)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/d217qr28g)**||
|| Central African Republic: |[2001-2021](https://scholar.colorado.edu/concern/datasets/pv63g1576) |**[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/zp38wf351)**||
|| Chad| [2001-2021](https://scholar.colorado.edu/concern/datasets/707958762)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/s4655j45w)**||
|| Comoros||**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/s7526f22t)**||
|| Democratic Republic of the Congo| [2001-2021](https://scholar.colorado.edu/concern/datasets/5425kb88g)|**[Nov 2000 - December 2024 (s1t1)](https://scholar.colorado.edu/concern/datasets/4b29b780c)**||
|| Djibouti| [2001-2021](https://scholar.colorado.edu/concern/datasets/1831cm01x)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/rf55z935n)**||
|| Egypt||**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/vm40xt27j)**||
|| Equatorial Guinea| [2001-2021](https://scholar.colorado.edu/concern/datasets/vx021g32b) | **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/j67315565)**||
|| Eritrea| [2001-2021](https://scholar.colorado.edu/concern/datasets/5m60qt182) | **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/jd472x93p)**||
|| eSwatini| [2001-2021](https://scholar.colorado.edu/concern/datasets/9w0324116)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2j62s642s)**||
|| Ethiopia| [2001-2021](https://scholar.colorado.edu/concern/datasets/z316q2977)| **[Nov 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/0c483m44k)**||
|| Gabon|[2001-2021](https://scholar.colorado.edu/concern/datasets/2z10wr67h)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5x21th08t)**||
|| The Gambia| [2001-2021](https://scholar.colorado.edu/concern/datasets/pn89d7911)|**[Nov 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/3t945s634)**||
|| Ghana|[2001-2021](https://scholar.colorado.edu/concern/datasets/2r36tz735)|**[Nov 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/79408001x)**||
|| Guinea|[2001-2021](https://scholar.colorado.edu/concern/datasets/05741s910)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/02870x464)**||
|| Guinea-Bissau| [2001-2021](https://scholar.colorado.edu/concern/datasets/nc580n858)|**[Nov 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/k643b3048)**||
|| Ivory Coast|[2001-2021](https://scholar.colorado.edu/concern/datasets/vq27zp62f)|**[Nov 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/rj430651d)**||
|| Kenya| [2001-2021](https://scholar.colorado.edu/concern/datasets/1j92g871c) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/v405sc03w)**||
|| Lesotho| [2001-2021](https://scholar.colorado.edu/concern/datasets/cr56n229w) | **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/w6634543v)**||
|| Liberia|[2001-2021](https://scholar.colorado.edu/concern/datasets/6h440t58k) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qv33rz453)**||
|| Libya||**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/z890rv992)**||
|| Madagascar| [2001-2021](https://scholar.colorado.edu/concern/datasets/fb494955x) |**[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/pv63g217f)**||
|| Malawi|[2001-2021](https://scholar.colorado.edu/concern/datasets/5999n464m) |**[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/t722hb67t)**||
|| Mali|[2001-2021](https://scholar.colorado.edu/concern/datasets/pr76f4544)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/td96k433q)**||
|| Mauritania|[2001-2021](https://scholar.colorado.edu/concern/datasets/x059c864s)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/v979v479r)**||
|| Mauritius|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1831cm54z)**||
|| Morocco| [2001-2021](https://scholar.colorado.edu/concern/datasets/td96k3751)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qj72p8679)**||
|| Mozambique|[2001-2021](https://scholar.colorado.edu/concern/datasets/1n79h5504)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2n49t338q)**||
|| Namibia|[2001-2021](https://scholar.colorado.edu/concern/datasets/db78td244)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5999n508g)**||
|| Niger| [2001-2021](https://scholar.colorado.edu/concern/datasets/m039k605q) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qn59q5683)**||
|| Nigeria|[2001-2021](https://scholar.colorado.edu/concern/datasets/cv43nx78p)|**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/j6731584b)**||
|| Republic of the Congo|[2001-2021](https://scholar.colorado.edu/concern/datasets/nk322f305)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/t435gf783)**||
|| Rwanda| [2001-2021](https://scholar.colorado.edu/concern/datasets/st74cr782)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/jw827d403)**||
|| Sao Tome and Principe||**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2227mr49g)**||
|| Senegal| [2001-2021](https://scholar.colorado.edu/concern/datasets/tt44pp176) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/v118rg081)**||
|| Sierra Leone|[2001-2021](https://scholar.colorado.edu/concern/datasets/5712m779r) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/r781wh64b)**||
|| Somalia| [2001-2021](https://scholar.colorado.edu/concern/datasets/xd07gt798) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/q524jq68k)**||
|| Somaliland| [2001-2021](https://scholar.colorado.edu/concern/datasets/8c97kr53f) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/08612q175)**||
|| South Africa| [2001-2021](https://scholar.colorado.edu/concern/datasets/rf55z8833) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/z603r0092)**||
|| South Sudan| [2001-2021](https://scholar.colorado.edu/concern/datasets/b2773w89g) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1g05fd479)**||
|| Sudan| [2001-2021](https://scholar.colorado.edu/concern/datasets/g158bj37v) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/9880vs542)**||
|| Tanzania| [2001-2021](https://scholar.colorado.edu/concern/datasets/7w62f947x) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/76537320w)**||
|| Togo |[2001-2021](https://scholar.colorado.edu/concern/datasets/fj236325p) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qf85nd104)**||
|| Tunisia|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/kp78gj27x)**||
|| Uganda| [2001-2021](https://scholar.colorado.edu/concern/datasets/hh63sx004)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/7d278v62s)**||
|| Zambia| [2001-2021](https://scholar.colorado.edu/concern/datasets/6108vc441) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qv33rz46c)**||
|| Zimbabwe| [2001-2021](https://scholar.colorado.edu/concern/datasets/f7623d95c)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/np193b769)**||
| Central North Africa |Libya, Algeria, Tunisia| [2001-2021](https://scholar.colorado.edu/concern/datasets/8910jv77j)|||
| **Asia** |||||
|| Afghanistan || **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/rj4306237)** ||
|| Armenia ||**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/7h149r64k)**||
|| Azerbaijan|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qb98mg90m)**||
|| Bhutan| [2001-2021](https://scholar.colorado.edu/concern/datasets/n009w342h) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2b88qd74k)**||
|| Bangladesh|[2001-2021](https://scholar.colorado.edu/concern/datasets/d791sh33k) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5d86p199z)**||
|| Bahrain|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/c534fq70v)**||
|| Cambodia|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2r36v0354)**||
|| China| [2001-2021](https://scholar.colorado.edu/concern/datasets/qz20st810) | **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/tt44pp77n)**||
|| India| [2001-2021](https://scholar.colorado.edu/concern/datasets/ht24wk47t)|**[Nov 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/cv43nz55k)**||
|| Israel|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5m60qt777)**||
|| Iraq|| **[Nov 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/w95052166)**||
|| Iran|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/00000172f)**||
|| Japan|[2001-2021](https://scholar.colorado.edu/concern/datasets/dz010r34v)|**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/zk51vj85q)**||
|| Kazakhstan|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1r66j289k)**||
|| Kuwait|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/z316q350j)**||
|| Kyrgyzstan|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/3197xn61c)**||
|| Lebanon|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qj72p8661)**
|| Laos| [2001-2021](https://scholar.colorado.edu/concern/datasets/bz60cx389)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/gq67js85b)**||
|| Macao|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/3f4627044)**||
|| Mongolia| [2001-2021](https://scholar.colorado.edu/concern/datasets/4x51hk21h) |**[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5t34sm20v)**||
|| Myanmar| [2001-2021](https://scholar.colorado.edu/concern/datasets/pk02cb86p)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/bv73c2197)**||
|| Nepal| [2001-2021](https://scholar.colorado.edu/concern/datasets/mk61rj10w)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/c534fq757)**||
|| North Korea| [2001-2021](https://scholar.colorado.edu/concern/datasets/3j333327g) |**[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/f7623f59q)**||
|| Oman|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/nz806134g)**||
|| Pakistan|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/8g84mn93n)**||
|| Palestine|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/nc580p37t)**||
|| Qatar|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/j098zc672)**||
|| Russia| [2001-2021](https://scholar.colorado.edu/concern/datasets/q811kk87t) |**[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/f7623f60g)**||
|| Saudi Arabia|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/g732db61d)**||
|| South Korea| [2001-2021](https://scholar.colorado.edu/concern/datasets/pg15bg177)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/kw52j990z)**||
|| Sri Lanka| [2001-2021](https://scholar.colorado.edu/concern/datasets/9z9030982)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/ws859h30c)**||
|| Syria|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/765373194)**||
|| Taiwan| [2001-2021](https://scholar.colorado.edu/concern/datasets/df65v9276) |**[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/6t053h714)**||
|| Tajikistan|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/9k41zg39z)**||
|| Thailand| [2001-2021](https://scholar.colorado.edu/concern/datasets/xs55md39h) |**[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qj72p868k)**||
|| Turkmenistan|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/4m90dx24j)**||
|| United Arab Emirates|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5138jg718)**||
|| Uzbekistan|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/ff365704w)**||
|| Vietnam| [2001-2021](https://scholar.colorado.edu/concern/datasets/h702q7566)| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/rn301342p)**||
|| Yemen|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qb98mg91w)**||
| Caucasus |Armenia, Azerbaijan, Georgia| [2001-2021](https://scholar.colorado.edu/concern/datasets/gf06g385j)|||
| Central Asia |Turkmenistan, Kazakhstan, Uzbekistan, Kyrgystan, Tajikistan, Afghanistan, and Pakistan| [2001-2021](https://scholar.colorado.edu/concern/datasets/47429b07v)|||
| Middle East |Saudi Arabia, Qatar, Oman, Yemen, United Arab Emirates, Iraq, Jordan, Syria, Israel, Palestine, Lebanon, Egypt| [2001-2021](https://scholar.colorado.edu/concern/datasets/5d86p139h)|||
| **Australia**|||||
|| Whole Country || **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/x059c9152)** ||
|| **(state by state)** ||||
||Tasmania| [2001-2021](https://scholar.colorado.edu/concern/datasets/c534fq19w)|||
|| Victoria| [2001-2021](https://scholar.colorado.edu/concern/datasets/2r36tz74f)|||
 || New South Wales + Capital Territory| [2001-2021](https://scholar.colorado.edu/concern/datasets/37720d85c)|||
 || Queensland| [2001-2021](https://scholar.colorado.edu/concern/datasets/cr56n230n)|||
 || South Australia| [2001-2021](https://scholar.colorado.edu/concern/datasets/fn107015p)|||
 || Western Australia| [2001-2021](https://scholar.colorado.edu/concern/datasets/k35695559)|||
 || Northern Territory| [2001-2021](https://scholar.colorado.edu/concern/datasets/bn9997900)|||
| **Oceania** ||
||Philippines| [2001-2021](https://scholar.colorado.edu/concern/datasets/7d278v06f) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/b8515q26w)**||
 || Papua New Guinea| [2001-2021](https://scholar.colorado.edu/concern/datasets/3r074w183) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/bk128c732)**||
 || East Timor| [2001-2021](https://scholar.colorado.edu/concern/datasets/j098zc184) |**[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/c247dt86h)**||
 || New Caledonia|| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/8w32r765q)**||
 || New Zealand| [2001-2021](https://scholar.colorado.edu/concern/datasets/9g54xj88f)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/h989r495c)**||
 || Malaysia| [2001-2021](https://scholar.colorado.edu/concern/datasets/fq977w13f)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/3484zj46h)**||
 || Brunei| [2001-2021](https://scholar.colorado.edu/concern/datasets/mp48sd92p)| **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/0v8382100)**||
 || Indonesia| [2001-2021](https://scholar.colorado.edu/concern/datasets/p2676w918)|**[Nov 2000 - December 2024](https://scholar.colorado.edu/concern/datasets/3b591b55c)**||
 || Samoa|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/kw52j9917)**||
 || Singapore|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/12579t91m)**||
 || Solomon Islands|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/6t053h70v)**||
 || Vanuatu|| **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/z316q351t)**||

## Installation

There are two main ways to install firedpy. Method 1 is to install locally from source and method 2 is to run it out of a Docker container.

### Method 1 - Install from source:

#### 1.1 Package Installer for Python (pip):

- **Step #1**: Install the Geospatial Data Abstraction Library (GDAL) on your machine, if it is not installed already. The installation method depends on your operating system (OS) and package manager (for Linux and macOS). Below are a subset of installation options (see the GDAL Documentation page for official downloads: https://gdal.org/en/stable/download.html).
  
    **macOS (brew)**
    ```bash
    brew install gdal
    ```
  
    **Linux (Debian or Debian-derivatives such as Ubuntu):**
    ```bash
    sudo apt update
    sudo apt install gdal-bin libgdal-dev python3-gdal
    ```
  
    **Linux (Fedora, RHEL, CentOS, Rocky Linux, and other OSs that use dnf)**
    ```bash
    sudo dnf install gdal gdal-devel python3-gdal
    ```
  
    **Linux (Arch)**
    ```bash
    sudo pacman -Syu
    sudo pacman -S gdal python-gdal
    ```

    **Windows**

    The most commonly recommended approach for installing GDAL with Windows is to use the Conda method (described below), but you may also use the OSGeo4W network installer. Go to the OSGeo4W landing page and download and run the installer: https://trac.osgeo.org/osgeo4w/. More detailed installation instructions are available on that site.

- **Step 2**: Create a Python environment and activate it. There are many options for doing this, one way is to use Python's virtual environment package as follows:

    **macOS or Linux**:
    ```bash
    mkdir ~/envs
    python3 -m venv ~/envs/firedpy
    source ~/envs/firedpy/bin/activate
    ```
  
    **Windows (CMD)**:
    ```cmd
    mkdir ~/envs
    python3 -m venv ~/envs/firedpy
    ~/envs/firedpy/Scripts/activate.bat
    ```

- **Step #3**: Now you may install FIREDpy from source using pip. Clone this repository, change directories into it, and install using the pip installation command:
  
    ```bash
    git clone https://github.com/earthlab/firedpy.git
    cd firedpy
    pip install .
    ```

#### 1.2 Conda: 

  Anaconda and its package manager `conda` is a system-independent method and should work for macOS, Linux, or Windows.

  - **Step #1**: Install Anaconda
  The installation method you use will depend on your operating system. Visit Conda's installation site, choose either `Miniconda`, `Anaconda Distribution`, or `Miniforge` (all options should work well) and following the installation instructions there: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html. The easiest way to find an installer is through this miniconda repository: https://repo.anaconda.com/miniconda/. 

  - **Step #2**: Clone this repository to a local folder and change directories into it:

    ```bash
    git clone https://github.com/earthlab/firedpy.git
    cd firedpy
    ```

  - **Step 3**: Use the `environment.yml` file to create a conda environment with all dependencies, including the system-level GDAL libraries.

    ```bash
    conda env create -f environment.yml
    conda activate firedpy
    ```
    - *Recommendation: be explicit about the environment location*
    
    ```bash
    conda env create -f environment.yml --prefix /path/to/firedenv
    conda activate /path/to/firedenb
    ```

  - **Step 4**: Install FIREDpy with pip.

    ```bash
    pip install .
    ```


### Method 2 - Run from a Docker Container:

#### 2.1 Get the docker container running:

Note, the docker container has changed from `earthlab/firedpy` to `earthlabcu/firedpy`
  - Run the docker container in a detached state (-d) and bind it to an available port on localhost (-p 127.0.0.1:0:7681)
  - `docker run -d -p 127.0.0.1:0:7681 earthlabcu/firedpy:latest`
 
  - Call `docker ps` to get the name of the docker container you just created and the port it is running on.
  
  - Then get into the docker container by either running docker exec:

    `docker exec -it <silly_name> /bin/bash`

  - Or access the CLI from your browser. The output from docker ps will look like this:

|CONTAINER ID |  IMAGE |                     COMMAND                |   CREATED    |     STATUS       |         PORTS                |                         NAMES|
|-------------|----------------------------|------------------------|---------------|------------------|-----------------------------|------------------------------|
|58a8a6ed926a |  earthlabcu/firedpy:latest | "/bin/entry.sh ttyd …" | 2 minutes ago |  Up 2 minutes  | 127.0.0.1:32768->7681/tcp     |                stupefied_hypatia|


In this example the container is running on the host machine at 127.0.0.1:32768. It may be different when you run it. Access this location in your browser by copy and pasting it into your browser's address bar

#### 2.2 Copy firedpy outputs to your local machine

After creating a new fire product, it might be useful to get it out of the docker container in order to use it.

  - First, exit the docker container by typing

    `exit`

  - Second, copy the file out. Here we will use the example of a container with the name "unruffled_clarke". The `docker cp` command uses the syntax `docker cp <source> <destination>`. Files inside of a docker container will have a prefix of the docker container name (or container ID) followed by a colon, then with a normal path.

    Here is an example command using the container name:

    `docker cp unruffled_clarke:/home/firedpy/proj/outputs/shapefiles/fired_events_s5_t11_2020153.gpkg /home/Documents/fired_events_s5_t11_2020153.gpkg`

    Another example command using the container ID:

    `docker cp fa73c6d3e007:/home/firedpy/proj/outputs/shapefiles/fired_events_s5_t11_2020153.gpkg /home/Documents/fired_events_s5_t11_2020153.gpkg`


## Use:

  - To avoid having to type in your username and password every time, follow the instructions below to create a .netrc file in your home directory

`https://nsidc.org/data/user-resources/help-center/creating-netrc-file-earthdata-login`

  - Run the command-line interface `firedpy` or `firedpy --help` in your terminal for its help file:

  ```bash
  firedpy --help
  ```

  You should see the help file with a description of the CLI and each of its options:

  ```bash

  Usage: firedpy [OPTIONS]

    firedpy command-line interface.

  Options:
    --version                       Show the version and exit.
    -i, --interactive               Interactive Mode. Firedpy will prompt the user for all parameter
                                    values and confirm selections before proceeding. Defaults to
                                    False.
    -p, --project_directory TEXT    Interactive Mode. Firedpy will prompt the user for all parameter
                                    values and confirm selections before proceeding. Defaults to
                                    False.  [required]
    -n, --project_name TEXT         A name used to identify the output files of this project. Defaults
                                    to None, which will use the name of the parent run directory.
  ...
  ```

  Then you may either run the CLI with each non-default parameter you wish define. You must define the project directory to help avoid writing outputs in your current working directory. You must also define either a `--tiles`, `--shape_file` or `--country` option for your study area. Be careful, running with defaults and a fire-intensive study area can be a very resource intensive operation:

  ```bash
  firedpy --project_directory ~/scratch/firedpy_finland --country Finland
  ```

  You may also run firedpy with the `--interactive` or `-i` flag to be prompted with input questions for each option:

  ```bash
  firedpy -i
  ```

  In this case, you will see a series of prompts starting with

  ```bash
  Using interactive mode, please enter parameter values as prompted. Press enter to accept defaults if availabile (value in square brackets)

  Project directory path. Inputs and outputs will all be written here. Required (Use "." for the present directory).: 
  ```

  The following set of commands demonstrates how to run `firedpy` with increasingly detailed control over the run options:

  - Change the spatial and temporal parameters of the model run:

    ```bash
    firedpy -project_directory ~/scratch/firedpy_finland --country Finland --spatial_param 6 --temporal_param 10
    ```

  - Specify specific MODIS tiles instead of a full country study area:

    ```bash
    firedpy -p ~/scratch/firedpy_nw_brazil --tiles "h11v09 h12v09" -sp 6 -tp 10
    ```

  - Specify specific years to run (defaults to 2000 to 2025):

    ```bash
    firedpy -p ~/scratch/firedpy_nw_brazil -t "h11v09 h12v09" -sp 6 -tp 10 --start_year 2020 --end_year 2024
    ```

  - Add the most common level 3 Ecoregion as an attribute to each event:

    ```bash
    firedpy -p ~/scratch/firedpy_nw_brazil -t "h11v09 h12v09" -sp 6 -tp 10 -y1 2020 -y2 2024 --eco_region_type 3
    ```

  - Add landcover information and produce the daily burn file:

    ```bash
    firedpy -p ~/scratch/firedpy_nw_brazil -t "h11v09 h12v09" -sp 6 -tp 10 -y1 2020 -y2 2024 --et 3 --land_cover_type 4 --daily
    ```


### Parameter table
  
| parameter | value(s)| example | description|
|:--------------|:----------|:-----|:---------|
| --spatial_param, -sp | integer | -sp 5 | Pixel radius for moving window, defaults to 5|
| --temporal_param, -tp | integer | -tp 11 | Day radius for moving window, defaults to 11|
| --tiles, -t | character (MODIS tile) | -t h11v09 | which modis tiles should be used |
| --shape_file, -sf | character (shapefile) | -sf ~/firedpy/data/individual_countries/canada.gpkg |  Finds MODIS tiles to download based on the polygon and clips output by the shapefile boundaries
| --project_directory, -p | character| -p /home/firedpy/proj | Sets the output directory. Will also download all input data to this path |
| --eco_region_type, -et | character | -et na | Type of ecoregion, either 'world' or 'na'|
| --eco_region_level, -el | integer | -el 3 | If ecoregion type = na, the level (1-3) of North American ecoregions |
| --land_cover_type, -lt | integer and character | -lt 2 | Number (1-3) corresponding with a MODIS/Terra+Aqua Land Cover (MCD12Q1) category. You will need to also make an account at https://urs.earthdata.nasa.gov/home |
| --shape_type, -st | character | -st gpkg | Option to build a shapefile for the fired event as a GeoPackage (gpkg), ESRI shapefile (shp), or both  |
| --project_name, -n | character | -n fired_colorado | Specifies a base name for output files for the tables and shapefile outputs |
| --daily, -d | character (yes or no) | -d | Creates daily polygons, if not provided only event-level perimeters will be created |
| --start_year, -y1 |integer | -y1 2001 | Gets the HDF4 files from the MODIS tiles starting in this year. The first year avalible is 2001 |
| --end_year |integer | -y2 2021 | Gets the HDF4 files from the MODIS tiles ending in this year. The last year avalible is 2021 |
 
 
 
### Boundary files are available for use as areas of interest
 
 - Country boundaries are in **data/individual_countries**
 - Continent boundaries are in **data/continents**
 - United States state boundaries for the United States of America are in **dat/us_states**
 - Australian state boundaries are in **data/australian_states**
 - For example `firedpy -sf ~/firedpy/firedpy/data/us_states/colorado.gpkg` gets Colorado, USA. 
 - If using the ineractive option, when prompted for the name of the country, spaces will be replaced with spaces and case does not matter. 


## How to update the docker container

- step 0.1. install docker (go to the docker site for OS-specific instructions.)
- step 0.2. get a dockerhub account
- step 1. login to docker hub via the command-line
   - `docker login` or `sudo docker login`
- step 2. get the existing docker image set up
   - docker run -t -d earthlab/firedpy
- step 3. update from GitHub
   - git pull 
- step 4. build the docker container
   - `docker build -t earthlab/firedpy:latest .`
- step 5. **ENSURE THE SOFTWARE STILL WORKS BEFORE PUSHING**
   - `firedpy -aoi /home/firedpy/ref/individual_countries/samoa.gpkg`
- step 6. push it up to dockerhub
   - `docker push earthlab/firedpy:latest`
