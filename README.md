[![DOI](https://zenodo.org/badge/214283770.svg)](https://zenodo.org/badge/latestdoi/214283770)

# firedpy
A Python CLI for classifying fire events from the Collection 6 MODIS Burned Area Product.

This function will use a space-time window to classify individual burn detections into discrete events and output both a data table and shapefiles of these events. The user is able to specify the spatial and temporal parameters of the model, as well as the area of interest using either a shapefile or a list of MODIS Sinusoidal Projection tile IDs. Any area from the world may be selected, though in the current version memory constraints may limit the extent available for a single model run. This version is calibrated to handle the Conitiguous United States (CONUS) with a 16 GB machine, though work is underway to move more processing to disk for larger areas. Shapefiles include event-level and daily level polygons, providing both final and expanding event perimeters.


Note for windows users: create Anaconda environment from the provided environment file and ensure you can import gdal. The install firedpy from same directory. For example, assuming gdal is installed correctly:

### To Install:

  - Clone this repository into a local directory:
  
    `git clone https://github.com/earthlab/firedpy.git`
    `cd firedpy`

  - Establish and activate a conda environment:
  
    `conda env create -f environment.yaml`  
    `conda activate firedpy`  

  - Install locally:
  
    `python setup.py install`


### To Use:

  - Run firedpy with no options to download required data and a write data table of classified fire events to a temporary directory. This uses a CONUS as the default area of interest with a spatial parameter of 5 pixels (~2,316.6 km) and 11 days:
 
    `firedpy`
  
  - In your terminal use this command to print out the available options and descriptions:
  
    `firedpy --help`
  
  - Change the spatial and temporal parameters of the model run:
  
    `firedpy -spatial_param 6 -temporal_param 10`
   
  - Specify specific tiles and a local project_directory for required data and model outputs:
  
    `firedpy -spatial_param 6 -temporal_param 10 -tiles h11v09 h12v09 -proj_dir /home/user/fired_project`
    
  - Write event- and daily-level shapefiles as outputs in addition a data table:
  
    `firedpy -spatial_param 6 -temporal_param 10 -tiles h11v09 h12v09 -proj_dir /home/user/fired_project --shapefile`

