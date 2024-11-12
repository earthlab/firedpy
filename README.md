
[![DOI](https://zenodo.org/badge/214283770.svg)](https://zenodo.org/badge/latestdoi/214283770) [![Docker Automated build](https://img.shields.io/docker/automated/earthlab/firedpy?style=plastic)](https://hub.docker.com/repository/docker/earthlab/firedpy/builds) ![GitHub contributors](https://img.shields.io/github/contributors/earthlab/firedpy) [![GitHub issues](https://img.shields.io/github/issues/earthlab/firedpy)](https://github.com/earthlab/firedpy/issues) ![GitHub commit activity](https://img.shields.io/github/commit-activity/w/earthlab/firedpy) 

# FIREDpy - FIRe Event Delineation for python

A Python Command Line Interface for classifying fire events from the Collection 6 MODIS Burned Area Product.

This package uses a space-time window to classify individual burn detections from late 2001 to near-present into discrete events and return both a data table and shapefiles of these events. The user is able to specify the spatial and temporal parameters of the window, as well as the spatial and temporal extent, using either a shapefile or a list of MODIS Sinusoidal Projection tile IDs. Shapefiles include full event polygons by default, and the user has the option of having firedpy produce daily-level perimeters, providing a representation of both final and expanding event perimeters. 

Any area from the world may be selected. However, in the current version, memory constraints may limit the extent available for a single model run. Equatorial regions have much more fire activity, and may require much more RAM to process than a normal laptop will have.

More methodological information is at:

Balch, J. K., St. Denis, L. A., Mahood, A. L., Mietkiewicz, N. P., Williams, T. P., McGlinchy J,
and Cook, M. C. 2020. FIRED (Fire Events Delineation): An open, flexible algorithm & database
of U.S. fire events derived from the MODIS burned area product (2001-19). Remote
Sensing, 12(21), 3498; https://doi.org/10.3390/rs12213498

Description of the country-level data sets is at: 

Mahood, A.L. Lindrooth, E.J., Cook, M.C. and Balch, J.K. Country-level fire perimeter datasets (2001-2021). 2022. Nature Scientific Data, 9(458). https://doi.org/10.1038/s41597-022-01572-3

## Changes
 - 10/14/2024 FIREDpy V2.0
    - No longer using setup.py. See new instructions below for running it with Docker or installing it locally.
    - Improved fire grouping
    - Improved CLI
    - Access to MODIS burn area product Version 6.1 with support up to at least 2024 

### BUG ALERT: 

Many of the data products created in Fall 2021 may be shifted by a half pixel, and may lack a coordinate reference system. 

The problem is now fixed, so this will not affect new iterations of firedpy. We created a script, R/posthoc_fixes.R that contains a function to fix either or both of these problems.

Sometimes the server (fuoco.geog.umd.edu) that houses the MCD64A1 product used by firedpy is down. If this happens, you just need to wait until it comes back up.

See the issues tab for more bugs, or to report a new bug!

## Have you used firedpy?

The algorithm and derived data products are under active development. Please take this [survey](https://docs.google.com/forms/d/e/1FAIpQLSe7ycmS0HGIze2T6TIUox8hsu8nlGsxiUMww8SCeWHDZPhB-Q/viewform?usp=sf_link) to help us improve firedpy. Or just email admahood@gmail.com and Adam will be overjoyed to talk about firedpy.

## Current status of created products

Already-created products are linked below. They are housed in the CU Scholar data repository in the [Earth Lab Data collection](https://scholar.colorado.edu/collections/pz50gx05h), or [here](https://scholar.colorado.edu/catalog?f%5Bacademic_affiliation_sim%5D%5B%5D=Earth+Lab&locale=en). 

All of the created products have an event-level shapefile in .gpkg and .shp formats. Many countries also have the daily-level shapefile, but these were not created for most countries in Africa and Asia due to computational restrictions. 

## Click on a link below to download a fire perimeter dataset for your favorite country

### North America
 - Coterminous USA + Alaska: [2001-2021](https://scholar.colorado.edu/concern/datasets/d504rm74m) **[2001-2024](https://scholar.colorado.edu/concern/datasets/fx719p11c)**
 - US plus Canada: [2001-2021](https://scholar.colorado.edu/concern/datasets/8336h304x)
 - Canada: [2001-2021](https://scholar.colorado.edu/concern/datasets/gf06g388c) **[2001-2024](https://scholar.colorado.edu/concern/datasets/0p096866d)**
 - Hawaii: [2001-2021](https://scholar.colorado.edu/concern/datasets/7h149r06p)
 - Carribean (Barbados, Bahamas, Cayman Islands, Cuba, Dominican Republic, Haiti, Jamaica, Montserrat, Puerto Rico, Saint Kitts And Nevis, Trinidad And Tobago, British Virgin Islands, Guadeloupe, Saint Barthelemy): [2001-2021](https://scholar.colorado.edu/concern/datasets/x633f230f)
 - Mexico and Central America (Belize, Guatemala, Honduras, El Salvador, Nicaragua, Costa Rica, Panama): [2001-2021](https://scholar.colorado.edu/concern/datasets/vd66w1102)
 - Aruba: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/9880vs496)**
 - Montserrat: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/3r074w69k)**
 - Antigua and Barbuda: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/st74cs33x)**
 - The Bahamas **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/t722hb658)**
 - Barbados: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/m613n033c)**
 - Belize: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/q524jq679)**
 - British Virgin Islands: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/jq085m540)**
 - Cayman Islands: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/0r967557n)**
 - Cuba: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5712m8365)**
 - Curacao: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/ws859h28b)**
 - Costa Rica: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/mc87pr85g)**
 - Dominica: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/hd76s170r)**
 - Dominican Republic: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1j92g9289)**
 - El Salvador: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/m613n034n)**
 - Grenada: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/9880vs50z)**
 - Haiti: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/707959237)**
 - Honduras: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/n296x082p)**
 - Jamaica: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/xg94hr218)**
 - Mexico: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5x21th07j)**
 - Panama: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/6q182m92n)**
 - Puerto Rico: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/gt54kq07d)**
 - Saint Barthelemy: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/h989r4943)**
 - Saint Kitts and Nevis: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2f75r969z)**
 - Saint Martin: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/7d278v61h)**
 - Trinidad and Tobago: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/9c67wp39g)**
 - Turks and Caicos: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/cv43nz30q)**
 - US Virgin Islands: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/pv63g214m)**

### South America
 - Argentina: [2001-2021](https://scholar.colorado.edu/concern/datasets/5t34sk58k) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/td96k4315)**
 - Brazil: [2001-2021](https://scholar.colorado.edu/concern/datasets/05741s90q) **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/b8515q25m)**
 - Bolivia: [2001-2021](https://scholar.colorado.edu/concern/datasets/b2773w83t)
 - Chile: [2001-2021](https://scholar.colorado.edu/concern/datasets/qr46r2052) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/g732db604)**
 - Colombia: [2001-2021](https://scholar.colorado.edu/concern/datasets/mp48sd91d)
 - Ecuador: [2001-2021](https://scholar.colorado.edu/concern/datasets/pc289k34n) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5138jg700)**
 - Guyana: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5h73pz00k)**
 - Paraguay: [2001-2021](https://scholar.colorado.edu/concern/datasets/rb68xd05p)
 - Peru: [2001-2021](https://scholar.colorado.edu/concern/datasets/x346d5616) **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/wh246t78k)**
 - Suriname: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/n583xw60q)**
 - Uruguay: [2001-2021](https://scholar.colorado.edu/concern/datasets/q524jq130) **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/q237ht604)**
 - Venezuela: [2001-2021](https://scholar.colorado.edu/concern/datasets/7m01bm95m)
 - Northern South America (Suriname, French Guiana, Guyana): [2001-2021](https://scholar.colorado.edu/concern/datasets/qv33rx839)

 
[Entire Western hemisphere from Jan 2017 to March 2020, intended for use in conjunction with GOES16 active fire detections.](https://scholar.colorado.edu/concern/datasets/d217qq78g)

### Europe
 - Aland: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5138jg697)**
 - Albania: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/wp988m53f)**
 - Andorra: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/xp68kh87x)**
 - Belarus: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/r781wh632)**
 - Belgium: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/hh63sx489)**
 - Bosnia and Herzegovina **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/dz010s09q)**
 - Croatia: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/kh04dr144)**
 - Cyprus: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/tx31qk48f)**
 - Czechia: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/fq977w613)**
 - Estonia: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/73666649m)**
 - Finland: [2001-2021](https://scholar.colorado.edu/concern/datasets/6395w836j) **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/pk02cc33v)**
 - Germany: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/6w924d770)**
 - Greece: [2001-2021](https://scholar.colorado.edu/concern/datasets/bc386k355) **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/z316q349s)**
 - Greenland: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/bv73c218z)**
 - Guernsey: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/bz60cx90j)**
 - Italy: [2001-2021](https://scholar.colorado.edu/concern/datasets/v979v416g) **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/tb09j7441)**
 - Kosovo: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/f1881n864)**
 - Latvia: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/fn1070711)**
 - Liechtenstein: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/f7623f56w)**
 - Lithuania: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/qb98mg925)**
 - Macedonia: **[Nov 2000 - July  2024](https://scholar.colorado.edu/concern/datasets/g158bj86s)**
 - Malta: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/sx61dp12d)**
 - Moldova: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/d217qr276)**
 - Monaco: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2514nn122)**
 - Montenegro: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1g05fd436)**
 - Netherlands: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/0k225c86c)**
 - Northern Cyprus: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/4m90dx238)**
 - Norway: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/70795924h)**
 - Poland: [2001-2021](https://scholar.colorado.edu/concern/datasets/6969z264g) **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/4b29b7569)**
 - Portugal: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/f4752j497)**
 - Russia: [2001-2021](https://scholar.colorado.edu/concern/datasets/q811kk87t)
 - San Marino: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/zk51vj568)**
 - Serbia: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5h73pz01v)**
 - Slovakia: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/cr56n2793)**
 - Slovenia: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/tq57ns75v)**
 - Spain: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/xw42n953k)**
 - Sweden: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1j92g930b)**
 - Switzerland: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1g05fd42x)**
 - United Kingdom: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/sf268684c)**
 - UK and Ireland: [2001-2021](https://scholar.colorado.edu/concern/datasets/pc289k33c)
 - Spain & Portugal: [2001-2021](https://scholar.colorado.edu/concern/datasets/gb19f7006)
 - Western Europe (France, Germany, Poland, Switzerland, Belgium, Netherlands, Luxembourg and Austria): [2001-2021](https://scholar.colorado.edu/concern/datasets/v692t736f)
 - Central to Southern Europe (Estonia, Latvia, Lithuania, Belarus, Ukraine, Czech Republic, Slovakia, Hungary, Romania, Bulgaria, Montenegro, Bosnia, Turkey, Republic Of Moldova, Serbia, Albania, Slovenia, and North Macedonia): [2001-2021](https://scholar.colorado.edu/concern/datasets/7h149r07z)
 - Northern Europe (Iceland, Sweden, Norway, and Denmark): [2001-2021](https://scholar.colorado.edu/concern/datasets/sb397945f)


### Africa

 - Algeria: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/bc386k800)**
 - Angola: [2001-2021](https://scholar.colorado.edu/concern/datasets/t435gf21z)
 - Benin: [2001-2021](https://scholar.colorado.edu/concern/datasets/z603qz58m) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/fb494b180)**
 - Botswana: [2001-2021](https://scholar.colorado.edu/concern/datasets/b8515p69g) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/nv9354593)**
 - Burundi: [2001-2021](https://scholar.colorado.edu/concern/datasets/3f462659h) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/vt150k986)**
 - Burkina Faso: [2001-2021](https://scholar.colorado.edu/concern/datasets/9g54xj875) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/x633f2869)**
 - Cabo Verde: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/v692t782j)**
 - Cameroon: [2001-2021](https://scholar.colorado.edu/concern/datasets/x920fz208) **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/d217qr28g)**
 - Central North Africa (Libya, Algeria, Tunisia): [2001-2021](https://scholar.colorado.edu/concern/datasets/8910jv77j)
 - Chad: [2001-2021](https://scholar.colorado.edu/concern/datasets/707958762)
 - Central African Republic: [2001-2021](https://scholar.colorado.edu/concern/datasets/pv63g1576)
 - Comoros: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/s7526f22t)**
 - Democratic Republic of the Congo: [2001-2021](https://scholar.colorado.edu/concern/datasets/5425kb88g)
 - Djibouti: [2001-2021](https://scholar.colorado.edu/concern/datasets/1831cm01x) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/rf55z935n)**
 - Egypt:  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/vm40xt27j)**
 - Equatorial Guinea: [2001-2021](https://scholar.colorado.edu/concern/datasets/vx021g32b)  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/j67315565)**
 - Eritrea: [2001-2021](https://scholar.colorado.edu/concern/datasets/5m60qt182)  **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/jd472x93p)**
 - eSwatini: [2001-2021](https://scholar.colorado.edu/concern/datasets/9w0324116)
 - Ethiopia: [2001-2021](https://scholar.colorado.edu/concern/datasets/z316q2977)
 - Gabon: [2001-2021](https://scholar.colorado.edu/concern/datasets/2z10wr67h)
 - The Gambia: [2001-2021](https://scholar.colorado.edu/concern/datasets/pn89d7911)
 - Ghana: [2001-2021](https://scholar.colorado.edu/concern/datasets/2r36tz735)
 - Guinea: [2001-2021](https://scholar.colorado.edu/concern/datasets/05741s910)
 - Guinea-Bissau: [2001-2021](https://scholar.colorado.edu/concern/datasets/nc580n858)
 - Ivory Coast: [2001-2021](https://scholar.colorado.edu/concern/datasets/vq27zp62f)
 - Kenya: [2001-2021](https://scholar.colorado.edu/concern/datasets/1j92g871c)
 - Lesotho: [2001-2021](https://scholar.colorado.edu/concern/datasets/cr56n229w)
 - Liberia: [2001-2021](https://scholar.colorado.edu/concern/datasets/6h440t58k)
 - Libya: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/z890rv992)**
 - Madagascar: [2001-2021](https://scholar.colorado.edu/concern/datasets/fb494955x)
 - Malawi: [2001-2021](https://scholar.colorado.edu/concern/datasets/5999n464m)
 - Mali: [2001-2021](https://scholar.colorado.edu/concern/datasets/pr76f4544) **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/td96k433q)**
 - Mauritania: [2001-2021](https://scholar.colorado.edu/concern/datasets/x059c864s)
 - Mauritius: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/1831cm54z)**
 - Morocco: [2001-2021](https://scholar.colorado.edu/concern/datasets/td96k3751) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qj72p8679)**
 - Mozambique: [2001-2021](https://scholar.colorado.edu/concern/datasets/1n79h5504)
 - Namibia: [2001-2021](https://scholar.colorado.edu/concern/datasets/db78td244)
 - Niger: [2001-2021](https://scholar.colorado.edu/concern/datasets/m039k605q)
 - Nigeria: [2001-2021](https://scholar.colorado.edu/concern/datasets/cv43nx78p)
 - Republic of the Congo: [2001-2021](https://scholar.colorado.edu/concern/datasets/nk322f305)
 - Rwanda: [2001-2021](https://scholar.colorado.edu/concern/datasets/st74cr782)
 - Sao Tome and Principe: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2227mr49g)**
 - Senegal: [2001-2021](https://scholar.colorado.edu/concern/datasets/tt44pp176)
 - Sierra Leone: [2001-2021](https://scholar.colorado.edu/concern/datasets/5712m779r)
 - Somalia: [2001-2021](https://scholar.colorado.edu/concern/datasets/xd07gt798) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/q524jq68k)**
 - Somaliland: [2001-2021](https://scholar.colorado.edu/concern/datasets/8c97kr53f) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/08612q175)**
 - South Africa: [2001-2021](https://scholar.colorado.edu/concern/datasets/rf55z8833)
 - South Sudan: [2001-2021](https://scholar.colorado.edu/concern/datasets/b2773w89g)
 - Sudan: [2001-2021](https://scholar.colorado.edu/concern/datasets/g158bj37v)
 - Tanzania: [2001-2021](https://scholar.colorado.edu/concern/datasets/7w62f947x)
 - Togo: [2001-2021](https://scholar.colorado.edu/concern/datasets/fj236325p)
 - Tunisia: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/kp78gj27x)**
 - Uganda: [2001-2021](https://scholar.colorado.edu/concern/datasets/hh63sx004)
 - Zambia: [2001-2021](https://scholar.colorado.edu/concern/datasets/6108vc441)
 - Zimbabwe: [2001-2021](https://scholar.colorado.edu/concern/datasets/f7623d95c)

### Asia

 - Afghanistan: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/rj4306237)**
 - Armenia: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/7h149r64k)**
 - Azerbaijan: **[Nov 2000 - July 2024]( https://scholar.colorado.edu/concern/datasets/qb98mg90m)**
 - Bhutan: [2001-2021](https://scholar.colorado.edu/concern/datasets/n009w342h) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2b88qd74k)**
 - Bangladesh: [2001-2021](https://scholar.colorado.edu/concern/datasets/d791sh33k)
 - Bahrain: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/c534fq70v)**
 - Cambodia: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/2r36v0354)**
 - China: [2001-2021](https://scholar.colorado.edu/concern/datasets/qz20st810)
 - India: [2001-2021](https://scholar.colorado.edu/concern/datasets/ht24wk47t)
 - Israel: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5m60qt777)**
 - Japan: [2001-2021](https://scholar.colorado.edu/concern/datasets/dz010r34v)
 - Kuwait: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/z316q350j)**
 - Kyrgyzstan: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/3197xn61c)**
 - Lebanon: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qj72p8661)**
 - Laos: [2001-2021](https://scholar.colorado.edu/concern/datasets/bz60cx389)
 - Macao: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/3f4627044)**
 - Mongolia: [2001-2021](https://scholar.colorado.edu/concern/datasets/4x51hk21h)
 - Myanmar: [2001-2021](https://scholar.colorado.edu/concern/datasets/pk02cb86p)
 - Nepal: [2001-2021](https://scholar.colorado.edu/concern/datasets/mk61rj10w)
 - North Korea: [2001-2021](https://scholar.colorado.edu/concern/datasets/3j333327g) **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/f7623f59q)**
 - Oman: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/nz806134g)**
 - Palestine: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/nc580p37t)**
 - Qatar: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/j098zc672)**
 - Saudi Arabia: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/g732db61d)**
 - South Korea: [2001-2021](https://scholar.colorado.edu/concern/datasets/pg15bg177) **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/kw52j990z)**
 - Sri Lanka: [2001-2021](https://scholar.colorado.edu/concern/datasets/9z9030982) **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/ws859h30c)**
 - Taiwan: [2001-2021](https://scholar.colorado.edu/concern/datasets/df65v9276) **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/6t053h714)**
 - Tajikistan: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/9k41zg39z)**
 - Thailand: [2001-2021](https://scholar.colorado.edu/concern/datasets/xs55md39h)
 - Turkmenistan: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/4m90dx24j)**
 - United Arab Emirates: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/5138jg718)**
 - Uzbekistan: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/ff365704w)**
 - Vietnam: [2001-2021](https://scholar.colorado.edu/concern/datasets/h702q7566)
 - Yemen: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/qb98mg91w)**
 - Caucasus (Armenia, Azerbaijan, Georgia): [2001-2021](https://scholar.colorado.edu/concern/datasets/gf06g385j)
 - Central Asia (Turkmenistan, Kazakhstan, Uzbekistan, Kyrgystan, Tajikistan, Afghanistan, and Pakistan): [2001-2021](https://scholar.colorado.edu/concern/datasets/47429b07v)
 - Middle East (Saudi Arabia, Qatar, Oman, Yemen, United Arab Emirates, Iraq, Jordan, Syria, Israel, Palestine, Lebanon, Egypt): [2001-2021](https://scholar.colorado.edu/concern/datasets/5d86p139h)

### Australia

 - Whole Country: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/x059c9152)**

### Australia (state by state)

 - Tasmania: [2001-2021](https://scholar.colorado.edu/concern/datasets/c534fq19w)
 - Victoria: [2001-2021](https://scholar.colorado.edu/concern/datasets/2r36tz74f)
 - New South Wales + Capital Territory: [2001-2021](https://scholar.colorado.edu/concern/datasets/37720d85c)
 - Queensland: [2001-2021](https://scholar.colorado.edu/concern/datasets/cr56n230n)
 - South Australia: [2001-2021](https://scholar.colorado.edu/concern/datasets/fn107015p)
 - Western Australia: [2001-2021](https://scholar.colorado.edu/concern/datasets/k35695559)
 - Northern Territory: [2001-2021](https://scholar.colorado.edu/concern/datasets/bn9997900)

### Oceania

 - Philippines: [2001-2021](https://scholar.colorado.edu/concern/datasets/7d278v06f) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/b8515q26w)**
 - Papua New Guinea: [2001-2021](https://scholar.colorado.edu/concern/datasets/3r074w183)
 - East Timor: [2001-2021](https://scholar.colorado.edu/concern/datasets/j098zc184) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/c247dt86h)**
 - New Caledonia: **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/8w32r765q)**
 - New Zealand: [2001-2021](https://scholar.colorado.edu/concern/datasets/9g54xj88f) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/h989r495c)**
 - Malaysia: [2001-2021](https://scholar.colorado.edu/concern/datasets/fq977w13f) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/3484zj46h)**
 - Brunei: [2001-2021](https://scholar.colorado.edu/concern/datasets/mp48sd92p) **[Nov 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/0v8382100)**
 - Indonesia: [2001-2021](https://scholar.colorado.edu/concern/datasets/p2676w918)
 - Samoa: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/kw52j9917)**
 - Singapore: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/12579t91m)**
 - Solomon Islands: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/6t053h70v)**
 - Vanuatu: **[November 2000 - July 2024](https://scholar.colorado.edu/concern/datasets/z316q351t)**


## Installation

There are two ways to install firedpy. Method one is to run it out of a docker container, Method 2 is to install locally.

### Method 1. Run from a Docker Container:

#### 1.1 Get the docker container running:

Note, the docker container has changed from `earthlab/firedpy` to `earthlabcu/firedpy`
  - Run the docker container in a detached state (-d) and bind it to an available port on localhost (-p 127.0.0.1:0:7681)
  - `docker run -d -p 127.0.0.1:0:7681 earthlabcu/firedpy:latest`
 
  - Call `docker ps` to get the name of the docker container you just created and the port it is running on.
  
  - Then get into the docker container by either running docker exec:

    `docker exec -it <silly_name> /bin/bash`

  - Or access the CLI from your browser. The output from docker ps will look like this:
    CONTAINER ID   IMAGE                       COMMAND                  CREATED         STATUS                 PORTS                                         NAMES
    58a8a6ed926a   earthlabcu/firedpy:latest   "/bin/entry.sh ttyd …"   2 minutes ago   Up 2 minutes           127.0.0.1:32768->7681/tcp                     stupefied_hypatia

    In this example the container is running on the host machine at 127.0.0.1:32768. It may be different when you run it. Access this location in your browser by copy and pasting it into your browser's address bar

#### 1.2 Copy firedpy outputs to your local machine

After creating a new fire product, it might be useful to get it out of the docker container in order to use it.

  - First, exit the docker container by typing

    `exit`

  - Second, copy the file out. Here we will use the example of a container with the name "unruffled_clarke". The `docker cp` command uses the syntax `docker cp <source> <destination>`. Files inside of a docker container will have a prefix of the docker container name (or container ID) followed by a colon, then with a normal path.

    Here is an example command using the container name:

    `docker cp unruffled_clarke:/home/firedpy/proj/outputs/shapefiles/fired_events_s5_t11_2020153.gpkg /home/Documents/fired_events_s5_t11_2020153.gpkg`

    Another example command using the container ID:

    `docker cp fa73c6d3e007:/home/firedpy/proj/outputs/shapefiles/fired_events_s5_t11_2020153.gpkg /home/Documents/fired_events_s5_t11_2020153.gpkg`


### Method 2. Local Installation Instructions:

  - Clone this repository to a local folder and change directories into it:

    `git clone https://github.com/earthlab/firedpy.git`

    `cd firedpy`
    
  - Ensure your anaconda setup has **conda-forge**, **channel_priority** set to **strict**, and **update your conda**.

    `conda update conda --yes`
    `conda config --add channels conda-forge`
    `conda config --set channel_priority strict`
    
  - You must have all packages listed in the environment.yaml installed using 'conda install -c conda-forge <package_name>'

  - Create and activate a conda environment:

    `conda env create -f environment.yml`

    `conda activate firedpy`

## Use:
  - Run firedpy with no options to be prompted with input questions for each option/attribute
     
    `python bin/firedpy.py` or if running from Docker container, simply `firedpy` 
    
  - Or use the following commands in your command line to specify the options/attributes you would like:   

  - In your terminal use this command to print out the available options and their descriptions:

    `python bin/firedpy.py --help`

  - Run firedpy with the default option to download required data and write a data table of classified fire events to a temporary directory. This uses CONUS as the default area of interest with a spatial parameter of 5 pixels (~2.3 km) and 11 days:

    `python bin/firedpy.py --default`

  - Change the spatial and temporal parameters of the model run:

    `python bin/firedpy.py -spatial 6 -temporal 10`

  - Specify specific tiles and a local project_directory for required data and model outputs:

    `python bin/firedpy.py -spatial 6 -temporal 10 -aoi h11v09 h12v09 -proj_dir /home/<user>/fired_project`

  - Write shapefiles as outputs in addition to the data table:

    `python bin/firedpy.py -spatial 6 -temporal 10 -aoi h11v09 h12v09 -proj_dir /home/<user>/fired_project --shapefile`

  - Add the most common level 3 Ecoregion as an attribute to each event:

    `python bin/firedpy.py bin/firedpy.py -spatial 6 -temporal 10 -aoi h11v09 h12v09 -proj_dir /home/<user>/fired_project --shapefile -ecoregion_level 3`

  - Add landcover information and produce the daily burn file

    `python bin/firedpy.py -spatial 6 -temporal 10 -aoi h11v09 h12v09 -proj_dir /home/<user>/fired_project --shapefile -ecoregion_level 3 -landcover_type 1 -daily yes`

  For more information about each parameter, use:

    'python bin/firedpy.py --help'
    
    
### Parameter table (under construction)
  
| parameter | value(s)| example | description|
|:--------------|:----------|:-----|:---------|
| -spatial | integer | -spatial 5 | pixel radius for moving window, defaults to 5|
| -temporal | integer | -temporal 11 | day radius for moving window, defaults to 11|
| -aoi | character (MODIS tile) | -aoi h11v09 | which modis tiles should be used |
| -aoi | character (shapefile) | -aoi /home/firedpy/individual_countries/canada.gpkg | figures out which modis tiles to download based on the polygon -- **polygon must be in the same projection as MODIS MCD64** -- all the polygons in the *ref* folder are correctly projected and can be used as crs templates to prepare other polygons. |
| -proj_dir| character| -proj_dir /home/firedpy/proj | which directory should firedpy operate within? Defaults to a folder called "proj" within the current working directory.|
| -ecoregion_type | character | -ecoregion_type na | type of ecoregion, either world or na|
 | -ecoregion_level | integer | -ecoregion_level 3 | if ecoregion type = na, the level (1-3) of North American ecoregions |
 | -landcover_type | integer and character | -landcover_type 2:username:password | number (1-3) corresponding with a MODIS/Terra+Aqua Land Cover (MCD12Q1) category. You will need to also make an account at https://urs.earthdata.nasa.gov/home and include your login information within the argument. |
 | -shp_type | character | -shp_type gpkg | option to build a shapefile for the fired event in gpkg, ESRI shapefile (shp), both, or none  |
 | -file | character | -file fired_colorado | specifies the base of the file name for the tables and shapefile outputs, defaults to "fired", in the format: "(-file aruguement)_toYYYYDDD_(either events or daily).gpkg", with YYYY being the year, and DDD being the julian day of the last month in the time series. The example would output fired_colorado_to2021031_events.gpkg.|
 | -daily | character (yes or no) | -daily yes | creates daily polygons, if no just the event-level perimeters will be created. Defaults to no. |
 | -start_yr |integer | -start_yr 2001 | gets the hdf files from the MODIS tiles starting in this year. The first year avalible is 2001 |
 | -end_yr |integer | -end_yr 2021 | gets the hdf files from the MODIS tiles ending in this year. The last year avalible is 2021 |
 
 
 
### Boundary files are available for use as areas of interest
 
 - Country boundaries are in **ref/individual_countries**
 - Continent boundaries are in **ref/continents**
 - United States state boundaries for the United States of America are in **ref/us_states**
 - Australian state boundaries are in **ref/australian_states**
 - For example `python bin/firedpy.py -aoi /home/firedpy/ref/us_states/colorado.gpkg`, and so on. Every space is a '_'. 
 - If using the user input option, when prompted for the name of the continent, country, or state use "_" for spaces. 
 - **Ensure that the input shapefiles are in the modis sinusiodal projection**


## How to update the docker container

- step 0.1. install docker (go to the docker website for OS-specific instructions.)
- step 0.2. get a dockerhub account
- step 1. login to docker hub via the command line
   - `docker login` or `sudo docker login`
- step 2. get the existing docker image set up
   - docker run -t -d earthlab/firedpy
- step 3. update from github
   - git pull 
- step 4. build the docker container
   - `docker build -t earthlab/firedpy:latest .`
- step 5. **ENSURE THE SOFTWARE STILL WORKS BEFORE PUSHING**
   - `firedpy -aoi /home/firedpy/ref/individual_countries/samoa.gpkg`
- step 6. push it up to dockerhub
   - `docker push earthlab/firedpy:latest`
