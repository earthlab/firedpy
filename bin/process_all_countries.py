import os
import subprocess as sp
import sys

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))

if __name__ == '__main__':
    countries = [c.strip('.gpkg') for c in os.listdir(os.path.join(PROJECT_DIR, 'ref', 'individual_countries'))]

    for country in countries:
        if not any([country in f for f in os.listdir(os.path.join(PROJECT_DIR, 'output', 'outputs'))]):
            process = sp.Popen([
                sys.executable,
                os.path.join(PROJECT_DIR, 'bin', 'firedpy.py'),
                "--full_csv", "n",
                "--n_cores", "1",
                "--tile_choice", "b",
                "--tile_name", country,
                "--daily", "y",
                "-spatial", "5",
                "-temporal", "11",
                "-shape_type", "both",
                "-land_cover_type", "1",
                "--cleanup", "y",
                "-start_year", "0",
                "-end_year", "0"
            ])  # Adjust shell=True if needed
            return_code = process.wait()  # Wait for process completion
