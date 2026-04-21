# -*- coding: utf-8 -*-
"""Run selected countries."""
import zipfile
from pathlib import Path

import pandas as pd

from firedpy import DATA_DIR
from firedpy.run import fired


PARAMETER_FPATH = DATA_DIR.joinpath("parameters_for_individual_countries.csv")
HOME = Path(__file__).parent.absolute()

# Override these for quick test runs; set TEST_COUNTRIES to None to process
# all countries, and PROJECT_NAME to None to fall back to using country names.
START_YEAR = 2000
END_YEAR = 2026
TEST_COUNTRIES = None  # e.g. ["spain", "portugal", "france"] or None
PROJECT_NAME = "v2025var"  # e.g. "v2025var"; if None, each country's name is used


def zip_run_outputs(directory, project_name, country, start_year, end_year,
                    spatial_param, temporal_param):
    """Zip the log, readme, and geofiles for a single firedpy run.

    Parameters
    ----------
    directory : pathlib.Path
        The project directory containing the ``outputs/`` subfolder.
    project_name : str
        The project name passed to :func:`firedpy.run.fired`.
    country : str
        Country name used as the area of interest label in the zip filename.
    start_year : int
        First year of the run.
    end_year : int
        Last year of the run.
    spatial_param : int | float
        Spatial parameter used for the run.
    temporal_param : int | float
        Temporal parameter used for the run.

    Returns
    -------
    pathlib.Path : Path to the created zip archive.
    """
    outputs_dir = directory / "outputs"
    run_name = f"fired_{project_name}"
    sp = int(spatial_param)
    tp = int(temporal_param)
    aoi = country.lower().replace(" ", "_")
    base = f"{run_name}_{aoi}_{start_year}-{end_year}_s{sp:02d}_t{tp:02d}"

    patterns = [
        f"{base}_*.log",
        f"{base}_readme.txt",
        f"{base}_events.*",
        f"{base}_daily.*",
    ]

    files_to_zip = []
    for pattern in patterns:
        files_to_zip.extend(outputs_dir.glob(pattern))

    zip_name = (
        f"fired_{project_name}_{aoi}_{start_year}-{end_year}"
        f"_s{sp:02d}_t{tp:02d}.zip"
    )
    zip_path = outputs_dir / zip_name

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for fpath in sorted(files_to_zip):
            zf.write(fpath, arcname=fpath.name)

    for fpath in files_to_zip:
        fpath.unlink()

    print(f"Zipped outputs to {zip_path}")
    return zip_path


def main():
    """Run firedpy for all countries."""
    # Table with all country names and fired parameters
    parameters = pd.read_csv(PARAMETER_FPATH)

    # To keep track of progress
    processed_fpath = HOME / "processed.txt"
    unprocessable_fpath = HOME / "unprocessable.txt"
    with open(unprocessable_fpath, 'r') as file:
        unprocessable = [f.strip('\n') for f in file.readlines()]
    with open(processed_fpath, 'r') as file:
        processed = [f.strip('\n') for f in file.readlines()]

    # Optional country filter for test runs
    if TEST_COUNTRIES is not None:
        parameters = parameters[parameters["country_name"].isin(TEST_COUNTRIES)]

    # Simple loop for now (reformat into slurm submission)
    for i, row in parameters.iterrows():
        country = row['country_name']
        sparam = row["spatial"]
        tparam = row["temporal"]
        run_name = PROJECT_NAME if PROJECT_NAME else country
        directory = HOME.joinpath("all")  # common directory so we don't download the same thing twice
        if country not in processed and country not in unprocessable:
            print(f"Processing {country}..")
            try:
                fired(
                    project_directory=directory,
                    project_name=run_name,
                    country=country,
                    start_year=START_YEAR,
                    end_year=END_YEAR,
                    spatial_param=row["spatial"],
                    temporal_param=row["temporal"],
                    daily=True,
                    shape_type="gpkg",
                    eco_region_type=None,
                    eco_region_level=None,
                    land_cover_type=1,
                    csv_type="none",
                    n_cores=0,
                    cleanup=False
                )
                with open(processed_fpath, "a") as file:
                    file.write(f"{country}\n")
                zip_run_outputs(
                    directory=directory,
                    project_name=run_name,
                    country=country,
                    start_year=START_YEAR,
                    end_year=END_YEAR,
                    spatial_param=sparam,
                    temporal_param=tparam,
                )
            except Exception as e:
                print(f"{country} failed: {e}")
                with open(unprocessable_fpath, 'a') as file:
                    file.write(f"{country}\n")


if __name__ == "__main__":
    main()
