import os
import subprocess as sp
import sys

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))

if __name__ == '__main__':
    countries = [c.strip('.gpkg') for c in os.listdir(os.path.join(PROJECT_DIR, 'ref', 'individual_countries'))]

    for country in countries:
        print(country)
        if country == "enya":
            country = "kenya"
        if country in ['mozambique', 'faroe_islands', 'western_sahara', 'french_polynesia']:
            continue
        if not any([country in f for f in os.listdir(os.path.join(PROJECT_DIR, 'output', 'outputs'))]):
            print("running", country)
            process = sp.Popen([
                sys.executable,
                os.path.join(PROJECT_DIR, 'bin', 'firedpy.py'),
                "--full_csv", "n",
                "--n_cores", "1",
                "--tile_choice", "b",
                "--tile_name", country,
                "--daily", "y",
                "-spatial", "1",
                "-temporal", "5",
                "-shape_type", "both",
                "-land_cover_type", "1",
                "--cleanup", "y",
                "-start_year", "0",
                "-end_year", "0"
            ])  # Adjust shell=True if needed
            return_code = process.wait()  # Wait for process completion
        else:
            print('Skipping', country)
