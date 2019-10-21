[![DOI](https://zenodo.org/badge/214283770.svg)](https://zenodo.org/badge/latestdoi/214283770)

# firedpy
A Python CLI for classifying fire events from the Collection 6 MODIS Burned Area Product.

Note for windows users: create Anaconda environment from the provided environment file and ensure you can import gdal. The install firedpy from same directory. For example, assuming gdal is installed correctly:

To Install:

  - Clone this repository into a local directory:
  
    `git clone https://github.com/earthlab/firedpy.git`
    `cd firedpy`

  - Establish and activate a conda environment:
  
    `conda env create -f environment.yaml`  
    `conda activate firedpy`  

  - Install locally:
  
    `python setup.py install`


To Use:

  - In your terminal use this command to print out the available options and descriptions:
  
    `firedpy --help`
  
  - Run firedpy with no options will download required data and model outputs for CONUS to a temporary directory:
    `firedpy`
  
