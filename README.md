[![DOI](https://zenodo.org/badge/214283770.svg)](https://zenodo.org/badge/latestdoi/214283770)
[![Build Status](https://travis-ci.com/earthlab/firedpy.svg?branch=master)](https://travis-ci.com/earthlab/firedpy)

# firedpy

A Python Command Line Interface for classifying fire events from the Collection 6 MODIS Burned Area Product.

This package uses a space-time window to classify individual burn detections from late 2001 to near-present into discrete events and return both a data table and shapefiles of these events. The user is able to specify the spatial and temporal parameters of the window, as well as the area of interest using either a shapefile or a list of MODIS Sinusoidal Projection tile IDs. Any area from the world may be selected. However, in the current version, memory constraints may limit the extent available for a single model run. This version is calibrated to handle the Conitiguous United States (CONUS) on a machine with 16 GB of RAM, though work is underway to move more processing to disk for larger areas. Shapefiles include full- and daily-level (still beta) event polygons, providing a representation of both final and expanding event perimeters.

There are two ways to install firedpy. Method 1 is to run out of a docker container, Method 2 is to install locally.

### Method 1. Run from a Docker Container:

#### 1.1 Get the docker container running:

    `docker run -t -d -p 8787:8787 earthlab/firedpy`
  
  - Call `docker ps` to get the name of the docker container you just created.
  
  - Then get into the docker container by running docker exec:

    `docker exec -it <silly_name> /bin/bash`
  
  - Then you will be inside of your docker container in the firedpy directory. Now, enter:
  
    `conda activate firedpy`
    
    And the environment should be ready to use.
    
#### 1.2 Copy firedpy outputs to your local machine

    After creating a new fire product, it might be useful to get it out of the docker container use it. 
    
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


### To Use:

  - In your terminal use this command to print out the available options and their descriptions:

    `firedpy --help`

  - Run firedpy with no options to download required data and write a data table of classified fire events to a temporary directory. This uses CONUS as the default area of interest with a spatial parameter of 5 pixels (~2.3 km) and 11 days:

    `firedpy`

  - Change the spatial and temporal parameters of the model run:

    `firedpy -spatial_param 6 -temporal_param 10`
 
  - Specify specific tiles and a local project_directory for required data and model outputs:

    `firedpy -spatial_param 6 -temporal_param 10 -tiles h11v09 h12v09 -proj_dir /home/<user>/fired_project`
  
  - Write shapefiles as outputs in addition to the data table:
  
    `firedpy -spatial_param 6 -temporal_param 10 -tiles h11v09 h12v09 -proj_dir /home/<user>/fired_project --shapefile`
  
  - Add the most common level 3 Ecoregion as an attribute to each event:
  
    `firedpy -spatial_param 6 -temporal_param 10 -tiles h11v09 h12v09 -proj_dir /home/<user>/fired_project --shapefile -ecoregion_level 3`
    
   Update 08/13/2020

