import os

from datetime import datetime, timedelta
from pathlib import Path

from firedpy import DATA_DIR


SUMMARY_TEMPLATE = DATA_DIR.joinpath("SUMMARY_TEMPLATE.txt")


def add_file_list(lines, output_directory):
    """Make a formatted list of files."""
    # Get all output files
    tables = list(output_directory.glob("**/*.csv"))
    gpkgs = list(output_directory.glob("**/*.gpkg"))
    shps = list(output_directory.glob("**/*.shp"))

    # Find the files line and add in the appropriate entries
    new_lines = []
    i = 1
    for line in lines:
        if "{files}" in line:
            if tables:
                new_lines.append(f"    1.{i} Tables\n")
                i += 1
                for item in tables:
                    new_lines.append(f"        - {str(item)}\n")
            if gpkgs:
                new_lines.append(f"    1.{i} Geopackages\n")
                i += 1
                for item in gpkgs:
                    new_lines.append(f"        - {str(item)}\n")
            if shps:
                new_lines.append(f"    1.{i} Shapefiles\n")
                i += 1
                for item in shps:
                    new_lines.append(f"        - {str(item)}\n")
        else:
            new_lines.append(line)

    return new_lines


def replace_values(parameters, line):
    """Replace instances of keys in line with their values.

    Parameters
    ----------
    paramters : dict
        A dictionary of parameters with keys representing strings to replace
        and values representing replacement strings.
    line : str
        A string that may contain the strings to be replaced.

    Returns
    -------
    str : A formatted string with replacements.
    """
    for string, replacement in parameters.items():
        line = line.replace(string, str(replacement))
    return line


def make_read_me(gdf, project_directory, tiles, spatial_param,
                 temporal_param, shapefile, runtime, n_cores,
                 peak_memory=None):
    """Write a summary file describing a firedpy run.

    Parameters
    ----------
    gdf : geopandas.geodataframe.GeoDataFrame
        The geodataframe output from a firedpy run.
    project_directory : str | pathlib.PosixPath
        The project directory containing all outputs of firedpy run used to
        build `gdf`.
    tiles : str
        The list of MODIS tile IDs used to build the fire event data.
    event_table_fpath : str
        The path to final event table.
    daily : bool
        If the firedpy run represents daily polygons or just the event-level
        perimeter for your analysis area.
    spatial_param : int
        The number of cells (~463 m resolution) to search for neighboring burn
        detections. Defaults to 5 cells in all directions.
    temporal_param : int
        The number of days to search for neighboring burn detections.
    shapefile : str | pathlib.PosixPath
        Path to a shapefile to use for the fire study area.
    runtime : int
        The total job runtime in minutes.
    n_cores : int
        The number of CPU cores used to perform the firedpy run.
    peak_memory : float
        The maximum memory usage reached during the firedpy run. Defaults to
        None, no summary written for this parameter.
    """
    # Infer the first and last date from the dataframe
    date1 = gdf["date"].min()
    date2 = gdf["date"].max()
    modis_date1 = date1.strftime("January %Y")
    modis_date2 = date2.strftime("December %Y")
    event_date1 = date1.strftime("%B %Y")
    event_date2 = date2.strftime("%B %Y")

    # List all output tables and spatial vector files
    project_directory = Path(project_directory).expanduser()
    output_directory = project_directory.joinpath("outputs")
    csvs = list(output_directory.glob("**/*csv"))

    # Assign specific files to variables

    # The run name is in the CSVs, ant there will always be at least one
    run_name = "_".join(csvs[0].name.split("_")[:-4])

    # Infer CPU count if default (0 for all) is used
    if n_cores == 0:
        n_cores = os.cpu_count() - 1

    # This is a stand in until we figure out how to track peak memory use
    if not peak_memory:
        peak_memory = "N/A"

    # Create a tile string
    tile_string = str(tiles)[1:-1]

    # Define replacement parameters
    parameters = {
        "{tile_string}": tile_string,
        "{modis_date1}": modis_date1,
        "{modis_date2}": modis_date2,
        "{spatial_param}": spatial_param,
        "{temporal_param}": temporal_param,
        "{run_name}": run_name,
        "{event_date1}": event_date1,
        "{event_date2}": event_date2,
        "{runtime}": f"{runtime:.2f}",
        "{peak_memory}": peak_memory,
        "{n_cores}": n_cores
    }

    # Read in the template
    with open(SUMMARY_TEMPLATE, "r") as template:
        lines = template.readlines()

    # Custom work for the files, which could be different each time
    lines_w_files = add_file_list(lines, output_directory)

    # Format template
    formatted_lines = []
    for line in lines_w_files:
        formatted_line = replace_values(parameters, line)
        formatted_lines.append(formatted_line)

    # Write to a README in the outputs directory
    fname = f"{run_name.upper()}_README.txt"
    fpath = output_directory.joinpath(fname)
    with open(fpath, "w") as summary:
        for line in formatted_lines:
            summary.write(line)
