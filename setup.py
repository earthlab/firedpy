# -*- coding: utf-8 -*-
import setuptools
#import subprocess as sp

# Issues with GDAL
#try:
#    gdal_version = sp.run(["gdalinfo",  "--version"],
#                          stdout=sp.PIPE)
#    gdal_version = gdal_version.stdout.decode("utf-8").replace("\n", "")
#    gdal_version = gdal_version[5:10]
#    pgvs = sp.run(["pip",  "install", "pygdal=="],
#                  stderr=sp.PIPE).stderr.decode("utf-8")
#    version_str = pgvs[pgvs.index("versions:") + 10: pgvs.index(")")]
#    versions = version_str.split(", ")
#    matching_vs = [v for v in versions if v[:5] == gdal_version]
#    pygdal_version = matching_vs[0]
#except ImportError:
#    print("Please install GDAL onto your system in order to use firedpy.")

# Description
with open("README.md", "r") as fh:
    long_description = fh.read()

print('find packages ', setuptools.find_packages())
setuptools.setup(
    name="firedpy",
    version="0.0.1",
    author="Travis Williams, Lise St. Denis",
    author_email="travis.williams@colorado.edu, lise.st.denis@colorado.edu",
    description="A CLI for classifying fires using a MODIS Burnt Area Product",
    license="MIT License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/earthlab/firedpy",
    platform=['any'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points = {'console_scripts': ['firedpy = firedpy.__main__:main']},
    packages=setuptools.find_packages(),
   install_requires = ['beautifulsoup4', 'bs4', 'dask', 'descartes', 'geopandas',
                       'lxml', 'netcdf4', 'numpy', 'pandas',
                       # 'pygdal==' + pygdal_version,
                       # 'pyyaml',
                       'scipy', 'toolz', 'tqdm', 'xarray',
                       'rasterio']
)
