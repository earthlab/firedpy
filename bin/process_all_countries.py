import os
import subprocess as sp
import sys

from bin.firedpy import cleanup_intermediate_files

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))

if __name__ == '__main__':
    countries = [c.strip('.gpkg') for c in os.listdir(os.path.join(PROJECT_DIR, 'ref', 'individual_countries'))]

    with open('unprocessable.txt', 'r') as file:
        # Append a new line to the file
        unprocessable = [f.strip('\n') for f in file.readlines()]
    with open('processed.txt', 'r') as file:
        # Append a new line to the file
        processed = [f.strip('\n') for f in file.readlines()]

    for country in countries:
        print(processed, unprocessable)
        if country not in processed and country not in unprocessable:
            print(country)
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
                "--cleanup", "n",
                "-start_year", "0",
                "-end_year", "0"
            ])  # Adjust shell=True if needed
            return_code = process.wait()  # Wait for process completion
            print(return_code)
            if return_code != 0:
                with open('unprocessable.txt', 'a') as file:
                    file.write(f"{country}\n")

                process = sp.Popen([
                    sys.executable,
                    'gocmd', 'mkdir', f'firedpy/{country}'
                ])
                return_code = process.wait()
                if return_code != 0:
                    continue
                process = sp.Popen([
                    sys.executable,
                    'gocmd', 'put', 'output/raster/burn_area/netcdfs', f'firedpy/{country}'
                ])  # Adjust shell=True if needed
                return_code = process.wait()
                if return_code != 0:
                    continue
                process = sp.Popen([
                    sys.executable,
                    'gocmd', 'put', 'output/raster/land_cover/mosaics', f'firedpy/{country}'
                ])  # Adjust shell=True if needed
                _ = process.wait()

                cleanup_intermediate_files(os.path.join(PROJECT_DIR, 'output'))
            else:
                with open('processed.txt', 'a') as file:
                    file.write(f"{country}\n")
