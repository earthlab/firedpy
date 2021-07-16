[![DOI](https://zenodo.org/badge/214283770.svg)](https://zenodo.org/badge/latestdoi/214283770)
[![Build Status](https://travis-ci.com/earthlab/firedpy.svg?branch=master)](https://travis-ci.com/earthlab/firedpy)

# firedpy - FIRe Event Delineation for python

A Python Command Line Interface for classifying fire events from the Collection 6 MODIS Burned Area Product.

This package uses a space-time window to classify individual burn detections from late 2001 to near-present into discrete events and return both a data table and shapefiles of these events. The user is able to specify the spatial and temporal parameters of the window, as well as the area of interest using either a shapefile or a list of MODIS Sinusoidal Projection tile IDs. Any area from the world may be selected. However, in the current version, memory constraints may limit the extent available for a single model run. This version is calibrated to handle the Contiguous United States (CONUS) on a machine with 16 GB of RAM, though work is underway to move more processing to disk for larger areas. Shapefiles include full- and daily-level event polygons, providing a representation of both final and expanding event perimeters.

More methodological information is at:

Balch, J. K., St. Denis, L. A., Mahood, A. L., Mietkiewicz, N. P., Williams, T. P., McGlinchy J,
and Cook, M. C. 2020. FIRED (Fire Events Delineation): An open, flexible algorithm & database
of U.S. fire events derived from the MODIS burned area product (2001-19). Remote
Sensing, 12(21), 3498; https://doi.org/10.3390/rs12213498

Already-created products:

#### Coterminous United States

 - events: https://scholar.colorado.edu/concern/datasets/nv935382p
 - daily: https://scholar.colorado.edu/concern/datasets/765372341

## Installation

There are two ways to install firedpy. Method one is to run it out of a docker container, Method 2 is to install locally.

### Method 1. Run from a Docker Container:

#### 1.1 Get the docker container running:

  - `docker run -t -d earthlab/firedpy`

  - Call `docker ps` to get the name of the docker container you just created.

  - Then get into the docker container by running docker exec:

    `docker exec -it <silly_name> /bin/bash`

  - Then you will be inside of the docker container in the firedpy directory. Now, enter:

    `conda activate firedpy`

    And the environment is ready to use.

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

  - Create and activate a conda environment:

    `conda env create -f environment.yaml`

    `conda activate firedpy`  

  - Install locally:

    `python setup.py install`


## Use:

  - In your terminal use this command to print out the available options and their descriptions:

    `firedpy --help`

  - Run firedpy with no options to download required data and write a data table of classified fire events to a temporary directory. This uses CONUS as the default area of interest with a spatial parameter of 5 pixels (~2.3 km) and 11 days:

    `firedpy`

  - Change the spatial and temporal parameters of the model run:

    `firedpy -spatial 6 -temporal 10`

  - Specify specific tiles and a local project_directory for required data and model outputs:

    `firedpy -spatial 6 -temporal 10 -tiles h11v09 h12v09 -proj_dir /home/<user>/fired_project`

  - Write shapefiles as outputs in addition to the data table:

    `firedpy -spatial 6 -temporal 10 -tiles h11v09 h12v09 -proj_dir /home/<user>/fired_project --shapefile`

  - Add the most common level 3 Ecoregion as an attribute to each event:

    `firedpy -spatial 6 -temporal 10 -tiles h11v09 h12v09 -proj_dir /home/<user>/fired_project --shapefile -ecoregion_level 3`

  - Add landcover information and produce the daily burn file

    `firedpy -spatial 6 -temporal 10 -tiles h11v09 h12v09 -proj_dir /home/<user>/fired_project --shapefile -ecoregion_level 3 -landcover_type 1 -daily yes`

  For more information about each parameter, use:

    'firedpy --help'
    
### Parameter table (under construction)
  
| parameter | value(s)| example | description|
|:--------------|:----------|:-----|:---------|
| -spatial | integer | -spatial 5 | pixel radius for moving window, defaults to 5|
| -temporal | integer | -temporal 11 | day radius for moving window, defaults to 11|
| -tiles | character | -tiles h11v09 | which modis tiles should be used|
|  |  | -tiles <polygon name> | figures out which modis tiles to download based on| the polygon |
| -proj_dir| character| -proj_dir /home/<user>/firedpy/proj | which directory should firedpy operate within? Defaults to the cwd.|
| -ecoregion_type | character | -ecoregion_type na | type of ecoregion, either world or na|
 | -ecoregion_level | integer | -ecoregion_level 3 | type of North American ecoregion |
 | -landcover_type | integer | -landcover_type 2 | number corresponding with a MODIS/Terra+Aqua Land Cover (MCD12Q1) category |
 | --shapefile | n/a | --shapefile | builds a .gpkg shapefile from the event data frame and saves it in the "outputs/shapefiles" folder in the project directory |
 | -file_name | character | -file_name fired_events.csv | saves the tables outputted as the designated file name, defaults to "fired" |
 | -daily | character | -daily yes | creates daily polygons, if no just the event-level perimeters will be created. Defaults to no. |
 
 
 
 ### Individual country boundaries available in ref/individual_countries
 
 Continent boundaries are in ref/continents
 
|Country Name                             |
|:----------------------------------------|
|afghanistan.gpkg                         |
|afghanistan.gpkg-shm                     |
|afghanistan.gpkg-wal                     |
|aland.gpkg                               |
|albania.gpkg                             |
|algeria.gpkg                             |
|american_samoa.gpkg                      |
|andorra.gpkg                             |
|angola.gpkg                              |
|anguilla.gpkg                            |
|antarctica.gpkg                          |
|antigua_and_barbuda.gpkg                 |
|argentina.gpkg                           |
|armenia.gpkg                             |
|aruba.gpkg                               |
|ashmore_and_cartier_islands.gpkg         |
|australia.gpkg                           |
|austria.gpkg                             |
|azerbaijan.gpkg                          |
|bahrain.gpkg                             |
|bangladesh.gpkg                          |
|barbados.gpkg                            |
|belarus.gpkg                             |
|belgium.gpkg                             |
|belize.gpkg                              |
|benin.gpkg                               |
|bermuda.gpkg                             |
|bhutan.gpkg                              |
|bolivia.gpkg                             |
|bosnia_and_herzegovina.gpkg              |
|botswana.gpkg                            |
|brazil.gpkg                              |
|british_indian_ocean_territory.gpkg      |
|british_virgin_islands.gpkg              |
|brunei.gpkg                              |
|bulgaria.gpkg                            |
|burkina_faso.gpkg                        |
|burundi.gpkg                             |
|cabo_verde.gpkg                          |
|cambodia.gpkg                            |
|cameroon.gpkg                            |
|canada.gpkg                              |
|cayman_islands.gpkg                      |
|central_african_republic.gpkg            |
|chad.gpkg                                |
|chile.gpkg                               |
|china.gpkg                               |
|colombia.gpkg                            |
|comoros.gpkg                             |
|cook_islands.gpkg                        |
|costa_rica.gpkg                          |
|croatia.gpkg                             |
|cuba.gpkg                                |
|curaçao.gpkg                             |
|cyprus.gpkg                              |
|czechia.gpkg                             |
|democratic_republic_of_the_congo.gpkg    |
|denmark.gpkg                             |
|djibouti.gpkg                            |
|dominica.gpkg                            |
|dominican_republic.gpkg                  |
|east_timor.gpkg                          |
|ecuador.gpkg                             |
|egypt.gpkg                               |
|el_salvador.gpkg                         |
|equatorial_guinea.gpkg                   |
|eritrea.gpkg                             |
|estonia.gpkg                             |
|eswatini.gpkg                            |
|ethiopia.gpkg                            |
|falkland_islands.gpkg                    |
|faroe_islands.gpkg                       |
|federated_states_of_micronesia.gpkg      |
|fiji.gpkg                                |
|finland.gpkg                             |
|france.gpkg                              |
|french_polynesia.gpkg                    |
|french_southern_and_antarctic_lands.gpkg |
|gabon.gpkg                               |
|gambia.gpkg                              |
|georgia.gpkg                             |
|germany.gpkg                             |
|ghana.gpkg                               |
|greece.gpkg                              |
|greenland.gpkg                           |
|grenada.gpkg                             |
|guam.gpkg                                |
|guatemala.gpkg                           |
|guernsey.gpkg                            |
|guinea-bissau.gpkg                       |
|guinea.gpkg                              |
|guyana.gpkg                              |
|haiti.gpkg                               |
|heard_island_and_mcdonald_islands.gpkg   |
|honduras.gpkg                            |
|hong_kong_s.a.r..gpkg                    |
|hungary.gpkg                             |
|iceland.gpkg                             |
|india.gpkg                               |
|indian_ocean_territories.gpkg            |
|indonesia.gpkg                           |
|iran.gpkg                                |
|iraq.gpkg                                |
|ireland.gpkg                             |
|isle_of_man.gpkg                         |
|israel.gpkg                              |
|italy.gpkg                               |
|ivory_coast.gpkg                         |
|jamaica.gpkg                             |
|japan.gpkg                               |
|jersey.gpkg                              |
|jordan.gpkg                              |
|kazakhstan.gpkg                          |
|kenya.gpkg                               |
|kiribati.gpkg                            |
|kosovo.gpkg                              |
|kuwait.gpkg                              |
|kyrgyzstan.gpkg                          |
|laos.gpkg                                |
|latvia.gpkg                              |
|lebanon.gpkg                             |
|lesotho.gpkg                             |
|liberia.gpkg                             |
|libya.gpkg                               |
|liechtenstein.gpkg                       |
|lithuania.gpkg                           |
|luxembourg.gpkg                          |
|macao_s.a.r.gpkg                         |
|macedonia.gpkg                           |
|madagascar.gpkg                          |
|malawi.gpkg                              |
|malaysia.gpkg                            |
|maldives.gpkg                            |
|mali.gpkg                                |
|malta.gpkg                               |
|marshall_islands.gpkg                    |
|mauritania.gpkg                          |
|mauritius.gpkg                           |
|mexico.gpkg                              |
|moldova.gpkg                             |
|monaco.gpkg                              |
|mongolia.gpkg                            |
|montenegro.gpkg                          |
|montserrat.gpkg                          |
|morocco.gpkg                             |
|mozambique.gpkg                          |
|myanmar.gpkg                             |
|namibia.gpkg                             |
|nauru.gpkg                               |
|nepal.gpkg                               |
|netherlands.gpkg                         |
|new_caledonia.gpkg                       |
|new_zealand.gpkg                         |
|nicaragua.gpkg                           |
|niger.gpkg                               |
|nigeria.gpkg                             |
|niue.gpkg                                |
|norfolk_island.gpkg                      |
|north_korea.gpkg                         |
|northern_cyprus.gpkg                     |
|northern_mariana_islands.gpkg            |
|norway.gpkg                              |
|oman.gpkg                                |
|pakistan.gpkg                            |
|palau.gpkg                               |
|palestine.gpkg                           |
|panama.gpkg                              |
|papua_new_guinea.gpkg                    |
|paraguay.gpkg                            |
|peru.gpkg                                |
|philippines.gpkg                         |
|pitcairn_islands.gpkg                    |
|poland.gpkg                              |
|portugal.gpkg                            |
|puerto_rico.gpkg                         |
|qatar.gpkg                               |
|republic_of_serbia.gpkg                  |
|republic_of_the_congo.gpkg               |
|romania.gpkg                             |
|russia.gpkg                              |
|rwanda.gpkg                              |
|saint_barthelemy.gpkg                    |
|saint_helena.gpkg                        |
|saint_kitts_and_nevis.gpkg               |
|saint_lucia.gpkg                         |
|saint_martin.gpkg                        |
|saint_pierre_and_miquelon.gpkg           |
|saint_vincent_and_the_grenadines.gpkg    |
|samoa.gpkg                               |
|san_marino.gpkg                          |
|são_tomé_and_principe.gpkg               |
|saudi_arabia.gpkg                        |
|senegal.gpkg                             |
|seychelles.gpkg                          |
|siachen_glacier.gpkg                     |
|sierra_leone.gpkg                        |
|singapore.gpkg                           |
|sint_maarten.gpkg                        |
|slovakia.gpkg                            |
|slovenia.gpkg                            |
|solomon_islands.gpkg                     |
|somalia.gpkg                             |
|somaliland.gpkg                          |
|south_africa.gpkg                        |
|south_georgia_and_the_islands.gpkg       |
|south_korea.gpkg                         |
|south_sudan.gpkg                         |
|spain.gpkg                               |
|sri_lanka.gpkg                           |
|sudan.gpkg                               |
|suriname.gpkg                            |
|sweden.gpkg                              |
|switzerland.gpkg                         |
|syria.gpkg                               |
|taiwan.gpkg                              |
|tajikistan.gpkg                          |
|tanzania.gpkg                            |
|thailand.gpkg                            |
|the_bahamas.gpkg                         |
|togo.gpkg                                |
|tonga.gpkg                               |
|trinidad_and_tobago.gpkg                 |
|tunisia.gpkg                             |
|turkey.gpkg                              |
|turkmenistan.gpkg                        |
|turks_and_caicos_islands.gpkg            |
|uganda.gpkg                              |
|ukraine.gpkg                             |
|united_arab_emirates.gpkg                |
|united_kingdom.gpkg                      |
|united_states_of_america.gpkg            |
|united_states_virgin_islands.gpkg        |
|uruguay.gpkg                             |
|uzbekistan.gpkg                          |
|vanuatu.gpkg                             |
|vatican.gpkg                             |
|venezuela.gpkg                           |
|vietnam.gpkg                             |
|wallis_and_futuna.gpkg                   |
|western_sahara.gpkg                      |
|yemen.gpkg                               |
|zambia.gpkg                              |
|zimbabwe.gpkg                            |
 

