# -*- coding: utf-8 -*-
import setuptools

# Description
with open("README.md", "r") as fh:
    long_description = fh.read()

print('find packages ', setuptools.find_packages())
setuptools.setup(
    name="firedpy",
    version="0.0.2",
    author="Travis Williams, Lise St. Denis, Adam Mahood, Maxwell C. Cook, Maxwell B. Joseph, Joseph McGlinchy, Estelle Lindrooth",
    author_email="adam.mahood@colorado.edu",
    description="A command line interface for classifying fires using a MODIS Burnt Area Product",
    license="MIT License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/earthlab/firedpy",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={'console_scripts': ['firedpy = firedpy.__main__:main']},
    packages=setuptools.find_packages()
)
