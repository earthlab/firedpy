# -*- coding: utf-8 -*-
"""
Fire Event Delineation for Python
"""
from pathlib import Path


PROJECT_DIR = Path(__file__).parent
DATA_DIR = PROJECT_DIR.joinpath("data")

__version__ = "2.0.0"

__doc__ = (
    """
    firedpy - Fire Event Delineation for Python
    -------------------------------------------
    A Python Command Line Interface for classifying fire events from the 
    Collection 6 MODIS Burned Area Product.
    """
)
