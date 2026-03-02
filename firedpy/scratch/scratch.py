# -*- coding: utf-8 -*-
"""Deprecated methods we're not willing to get rid of entirely yet.

Author: travis
Date: Fri Feb 13 06:08:46 PM MST 2026
"""
import os
import tarfile
import urllib.request

from glob import glob
from http.cookiejar import CookieJar


EARTHDATA_TEST_URL = (
    "https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.061/2019.01.01/"
    "BROWSE.MCD12Q1.A2019001.h10v09.061.2022169160720.1.jpg"
)


def test_earthdata_credentials(username, password):
    """Test access to `ers.earthdata.nasa.gov`

    Parameters
    ----------
    username : str
        Username for Earthdata account.
    password : str
        Password for Earthdata account.
    """
    # Earthdata Login
    password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    password_manager.add_password(None, "https://urs.earthdata.nasa.gov",
                                  username, password)

    # Create a cookie jar for storing cookies. This is used to store and return
    # the session cookie given to use by the data server (otherwise it will
    # just keep sending us back to Earthdata Login to authenticate).  Ideally,
    # we should use a file based cookie jar to preserve cookies between runs.
    # This will make it much more efficient.
    cookie_jar = CookieJar()

    # Install all the handlers
    opener = urllib.request.build_opener(
        urllib.request.HTTPBasicAuthHandler(password_manager),
        # urllib.request.HTTPHandler(debuglevel=1),  # Uncomment to see details
        # urllib.request.HTTPSHandler(debuglevel=1),  # of requests/responses
        urllib.request.HTTPCookieProcessor(cookie_jar))
    urllib.request.install_opener(opener)

    # Send a test URL
    request = urllib.request.Request(EARTHDATA_TEST_URL)
    urllib.request.urlopen(request)


def package_individual_country_outputs():
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
