[![DOI](https://zenodo.org/badge/214283770.svg)](https://zenodo.org/badge/latestdoi/214283770)
[![Build Status](https://travis-ci.com/earthlab/firedpy.svg?branch=master)](https://travis-ci.com/earthlab/firedpy)

# firedpy
A Python Command Line Interface for classifying fire events from the Collection 6 MODIS Burned Area Product.

This package will use a space-time window to classify individual burn detections from late 2001 to near-present into discrete events and return both a data table and shapefiles of these events. The user is able to specify the spatial and temporal parameters of the window, as well as the area of interest using either a shapefile or a list of MODIS Sinusoidal Projection tile IDs. Any area from the world may be selected. However, in the current version, memory constraints may limit the extent available for a single model run. This version is calibrated to handle the Conitiguous United States (CONUS) with a 16 GB machine, though work is underway to move more processing to disk for larger areas. Shapefiles include full- and daily-level event polygons, providing a representation of both final and expanding event perimeters.


### Installation instructions:

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

