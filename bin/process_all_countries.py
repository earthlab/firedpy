import os
import subprocess as sp
import sys

from bin.firedpy import cleanup_intermediate_files

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))

if __name__ == '__main__':
    countries = [c.rstrip('.gpkg') for c in os.listdir(os.path.join(PROJECT_DIR, 'ref', 'individual_countries'))]

    with open('unprocessable.txt', 'r') as file:
        # Append a new line to the file
        unprocessable = [f.strip('\n') for f in file.readlines()]
    with open('processed.txt', 'r') as file:
        # Append a new line to the file
        processed = [f.strip('\n') for f in file.readlines()]

    for country in countries:
        if country not in processed and country not in unprocessable:
            print(country)
            try:
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
                return_code = process.wait(timeout=8 * 60 * 60)  # Wait for process completion
                if return_code == 0:
                    with open('processed.txt', 'a') as file:
                        file.write(f"{country}\n")
                else:
                    with open('unprocessable.txt', 'a') as file:
                        file.write(f"{country}\n")
                    cleanup_intermediate_files(os.path.join(PROJECT_DIR, 'output'))
            except (sp.TimeoutExpired, MemoryError):
                # If the process exceeds the timeout, terminate it
                print("Process took too long, terminating...")
                process.kill()  # Forcefully terminate the process
                # stdout, stderr = process.communicate()
                # exit_code = process.returncode

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
