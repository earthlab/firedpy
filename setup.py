"""Install firedpy, requires GDAL."""
import subprocess as sp

from setuptools import setup


def get_gdal_version():
    """Return system GDAL version."""
    # This call should be available to Windows too, alternately `gdal-config`
    process = sp.Popen(
        ["gdalinfo", "--version"],
        stdout=sp.PIPE,
        stderr=sp.PIPE
    )
    sto, ste = process.communicate()
    if ste:
        raise OSError("GDAL is causing problems again. Make sure you can run "
                      "'gdalinfo --version' successfully in your terminal.")
    version = sto.decode().split()[1]
    return version


def get_requirements():
    """Get requirements and update gdal version number."""
    with open("requirements.txt", encoding="utf-8") as file:
        reqs = file.readlines()
    gdal_version = get_gdal_version()
    gdal_line = f"gdal=={gdal_version}.*\n"
    reqs.append(gdal_line)
    return reqs


setup(
    name="firedpy",
    packages=["firedpy"],
    install_requires=get_requirements(),
    include_package_data=True,
    package_data={
        "data": [
            "*"
        ]
    },
)
