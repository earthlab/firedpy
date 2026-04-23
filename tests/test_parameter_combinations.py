# -*- coding: utf-8 -*-
"""Test specific sites for basic functionality and particular issues.

Author: travis
Date: Sun Apr  5 09:32:37 AM MDT 2026
"""
import pytest

from pathlib import Path

from firedpy import DATA_DIR
from firedpy.run import fired


@pytest.mark.parametrize("cores", [0, 1])
@pytest.mark.parametrize("daily", [True, False])
@pytest.mark.parametrize("cleanup", [True, False])
def test_iceland(cores, daily, cleanup):
    """Test a sample of fires in Iceland for expected behavior."""
    # Run Fired for Iceland
    pdir = Path("./testruns_iceland").absolute()
    out = fired(
        project_directory=pdir,
        country="Iceland",
        project_name="testruns_iceland",
        start_year=2021,
        end_year=2022,
        spatial_param=8,
        temporal_param=3,
        daily=daily,
        shape_type="gpkg",
        eco_region_type="world",
        eco_region_level=3,
        land_cover_type=1,
        full_csv=True,
        n_cores=cores,
        cleanup=cleanup
    )

    # We expect three events in 2022
    n = out.shape[0]
    error_msg = f"Expected 3 burn detections for Iceland in 2022, got {n}."
    assert out.shape[0] == 3, error_msg


def test_shapefile():
    """Test a using a shapefile for the study area."""
    # Run Fired for Iceland
    pdir = Path("./testruns_shapefile").absolute()
    shp_fpath = DATA_DIR.joinpath("individual_countries/finland.gpkg")
    _ = fired(
        project_directory=pdir,
        country=None,
        project_name="testruns_shapefile",
        shape_file=shp_fpath,
        start_year=2000,
        end_year=2025,
        spatial_param=8,
        temporal_param=3,
        daily=True,
        shape_type="gpkg",
        eco_region_type=None,
        eco_region_level=3,
        land_cover_type=None,
        full_csv=True,
        n_cores=1,
        cleanup=True
    )
