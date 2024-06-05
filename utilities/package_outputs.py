import os
import tarfile
from glob import glob


def pacakge_individual_country_outputs():
    unique_outputs = list(set('_'.join(s.split('fired_')[1].split('.csv')[0].split('_')[:-1]) for s in
                              glob("output/fired*.csv")))

    shapefile_extensions = ['*.shp', '*.shx', '*.dbf', '*.prj', '*.sbn', '*.sbx', '*.fbn', '*.fbx', '*.ain', '*.aih',
                            '*.ixs', '*.mxs', '*.atx', '*.cpg']

    for unique_output in unique_outputs:
        csv_files = glob(os.path.join("output", f"fired_{unique_output}*.csv"))
        read_me = glob(os.path.join("output", "outputs", f"fired_{unique_output}*README.txt"))
        gpkg_files = glob(os.path.join("output", "outputs", "shapefiles", f"fired_{unique_output}*.gpkg"))
        shape_files = []
        for ext in shapefile_extensions:
            shape_files.extend(
                glob(os.path.join("output", "outputs", "shapefiles", f"fired_{unique_output}{ext}")))

        if csv_files and read_me and shape_files and gpkg_files:
            # Combine all file lists into one
            shape_file_package = csv_files + read_me + shape_files
            gpkg_file_package = csv_files + read_me + gpkg_files

            # Define the output tar file name
            shape_tar_file_name = os.path.join("output", f"{unique_output}_shp.tar")
            if os.path.exists(shape_tar_file_name):
                print(f'{shape_tar_file_name} already exists')
                continue
            with tarfile.open(shape_tar_file_name, "w") as tar:
                for file in shape_file_package:
                    tar.add(file, arcname=os.path.relpath(file, "output"))
            print(f"Created tar file: {shape_tar_file_name}")

            gpkg_tar_file_name = os.path.join("output", f"{unique_output}_gpkg.tar")
            if os.path.exists(gpkg_tar_file_name):
                print(f'{gpkg_tar_file_name} already exists')
                continue
            with tarfile.open(gpkg_tar_file_name, "w") as tar:
                for file in gpkg_file_package:
                    tar.add(file, arcname=os.path.basename(file))
            print(f"Created tar file: {gpkg_tar_file_name}")

            for file in csv_files + read_me + shape_files + gpkg_files:
                os.remove(file)
