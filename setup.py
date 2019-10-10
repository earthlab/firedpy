# -*- coding: utf-8 -*-
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="firedpy",
    version="0.0.1",
    author="Travis Williams, Lise St. Denis",
    author_email="travis.williams@colorado.edu, lise.st.denis@colorado.edu",
    description="A wildfire classification package",
    license="MIT License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/earthlab/firedpy",
    packages=setuptools.find_packages(),
    platform=['any'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
