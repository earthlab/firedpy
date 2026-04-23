# -*- coding: utf-8 -*-
"""Run selected countries."""
from pathlib import Path

import pandas as pd

from firedpy import DATA_DIR
from firedpy.run import fired


PARAMETER_FPATH = DATA_DIR.joinpath("parameters_for_individual_countries.csv")
HOME = Path(__file__).parent.absolute()
SLURM_TEMPLATE = (
"""#!/bin/bash

#SBATCH --account=PLACEHOLDER
#SBATCH --time=1:00:00
#SBATCH -o PLACEHOLDER.o
#SBATCH -e PLACEHOLDER.e
#SBATCH --job-name=<PLACEHOLDER>
#SBATCH --nodes=1
#SBATCH --mail-user=PLACEHOLDER
#SBATCH --mem=79000

echo Running on: $HOSTNAME, Machine Type: $MACHTYPE
echo CPU: $(cat /proc/cpuinfo | grep "model name" -m 1 | cut -d:  -f2)
echo RAM: $(free -h | grep  "Mem:" | cut -c15-21)

source ~/.bashrc
module load conda
conda activate /path/to/env/

python script.py
"""
)


def main():
    """Run firedpy for all countries."""
    # Table with all country names and fired parameters
    parameters = pd.read_csv(PARAMETER_FPATH)

    # To keep track of progress
    with open('unprocessable.txt', 'r') as file:
        unprocessable = [f.strip('\n') for f in file.readlines()]
    with open('processed.txt', 'r') as file:
        processed = [f.strip('\n') for f in file.readlines()]

    # Simple loop for now (reformat into slurm submission)
    for i, row in parameters.iterrows():
        country = row['country_name']
        sparam = row["spatial"]
        tparam = row["temporal"]
        directory = HOME.joinpath("all")  # have them all in a common directory so we don't download the same thing twice
        if country not in processed and country not in unprocessable:
            print(f"Processing {country}..")
            try:
                fired(
                    project_directory=directory,
                    project_name=country, # here we can add the s and t params to the filename
                    country=country,
                    start_year=2000,
                    end_year=2026,
                    spatial_param=row["spatial"],
                    temporal_param=row["temporal"],
                    daily=True,
                    shape_type="gpkg",
                    eco_region_type=None,
                    eco_region_level=None,
                    land_cover_type=1,
                    full_csv=False,
                    n_cores=0,
                    cleanup=False
                )
                with open("processed.txt", "a") as file:
                    file.write(f"{country}\n")
            except Exception as e:
                print(f"{country} failed: {e}")
                with open("unprocessable.txt", 'a') as file:
                    file.write(f"{country}\n")


if __name__ == "__main__":
    main()
