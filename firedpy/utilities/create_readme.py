import calendar
import logging
import os
import re

from pathlib import Path

from firedpy import DATA_DIR

logger = logging.getLogger(__name__)


def _parse_modis_date_range(nc_files):
    """Return (modis_date1, modis_date2) human-readable strings from NC paths.

    Each NC filename is expected to follow the month-aware naming convention:
    ``{tile}_{YYYY}-{MM}_{YYYY}-{MM}.nc`` (e.g. ``h18v02_2000-01_2026-03.nc``).

    Across all supplied NC files the function takes the earliest start and the
    latest end so that the reported range spans the full set of tiles.

    Returns a pair of formatted strings like ``"January 2000"`` and
    ``"March 2026"``, or ``None, None`` if no filenames match the pattern.
    """
    pat = re.compile(r"_(\d{4})-(\d{2})_(\d{4})-(\d{2})\.nc$")
    starts, ends = [], []
    for nc in nc_files:
        m = pat.search(str(nc))
        if m:
            sy, sm, ey, em = int(m.group(1)), int(m.group(2)), \
                             int(m.group(3)), int(m.group(4))
            starts.append((sy, sm))
            ends.append((ey, em))

    if not starts:
        return None, None

    min_sy, min_sm = min(starts)
    max_ey, max_em = max(ends)
    date1 = f"{calendar.month_name[min_sm]} {min_sy}"
    date2 = f"{calendar.month_name[max_em]} {max_ey}"
    return date1, date2


SUMMARY_TEMPLATE = DATA_DIR.joinpath("SUMMARY_TEMPLATE.txt")


def add_file_list(lines, output_directory, run_name, aoi, start_year,
                  end_year, spatial_param, temporal_param):
    """Make a formatted list of files produced by this specific run.

    Only includes files whose names match the run's naming pattern, so that
    re-running firedpy in the same output directory does not pollute the file
    list with outputs from previous runs.
    """
    sp = int(spatial_param)
    tp = int(temporal_param)
    aoi_part = f"_{aoi}" if aoi else ""
    run_prefix = (f"{run_name}{aoi_part}_{start_year}-{end_year}"
                  f"_s{sp:02d}_t{tp:02d}")

    tables = sorted(output_directory.glob(f"{run_prefix}*.csv"))
    gpkgs = sorted(output_directory.glob(f"{run_prefix}*.gpkg"))
    shps = sorted(output_directory.glob(f"{run_prefix}*.shp"))
    logs = sorted(output_directory.glob(f"{run_prefix}*.log"))
    readme_name = f"{run_prefix.lower()}_readme.txt"

    # Find the {files} placeholder and expand it
    new_lines = []
    i = 1
    for line in lines:
        if "{files}" in line:
            if tables:
                new_lines.append(f"    1.{i} Tables\n")
                i += 1
                for item in tables:
                    new_lines.append(f"        - {item.name}\n")
            if gpkgs:
                new_lines.append(f"    1.{i} Geopackages\n")
                i += 1
                for item in gpkgs:
                    new_lines.append(f"        - {item.name}\n")
            if shps:
                new_lines.append(f"    1.{i} Shapefiles\n")
                i += 1
                for item in shps:
                    new_lines.append(f"        - {item.name}\n")
            if logs:
                new_lines.append(f"    1.{i} Log files\n")
                i += 1
                for item in logs:
                    new_lines.append(f"        - {item.name}\n")
            new_lines.append(f"    1.{i} Run summary\n")
            new_lines.append(f"        - {readme_name}\n")
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
                 start_year, end_year, run_name=None, country=None,
                 peak_memory=None, aoi=None, nc_files=None):
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
    start_year : int
        The first year of fire events.
    end_year : int
        The last year of fire events.
    run_name : str | None
        The run name used to prefix output files (e.g. ``fired_v2.1``).
        Derived from the output directory name if not provided.
    country : str | None
        Country name used as the study area, if provided. Defaults to None.
    aoi: str | None
        Normalised area-of-interest label (e.g. country name, shapefile, tile),
        included in the readme filename.
    peak_memory : float
        The maximum memory usage reached during the firedpy run. Defaults to
        None, no summary written for this parameter.
    nc_files : list[str | pathlib.Path] | None
        Paths to the NetCDF cache files used in this run.  When supplied, the
        MODIS burned-area product date range is parsed from their filenames
        (``{tile}_{YYYY}-{MM}_{YYYY}-{MM}.nc``) and written to the abstract.
        Falls back to the fire-event date range when not provided or when the
        filenames do not match the expected pattern.
    """
    # Infer the first and last fire-event date from the dataframe
    date1 = gdf["date"].min()
    date2 = gdf["date"].max()
    event_date1 = date1.strftime("%B %Y")
    event_date2 = date2.strftime("%B %Y")

    # MODIS burned-area product date range from NC filenames
    modis_date1, modis_date2 = _parse_modis_date_range(nc_files or [])
    if modis_date1 is None:
        modis_date1, modis_date2 = event_date1, event_date2

    # Build a human-readable area-of-interest label:
    # custom shapefile stem > country name > tile IDs
    if shapefile:
        aoi_label = Path(shapefile).stem
    elif country:
        aoi_label = country
    else:
        aoi_label = str(tiles)[1:-1]

    # Set up output directory
    project_directory = Path(project_directory).expanduser()
    output_directory = project_directory.joinpath("outputs")

    # Use the provided run_name, or fall back to the project directory name
    if not run_name:
        run_name = project_directory.name

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
        "{event_date1}": event_date1,
        "{event_date2}": event_date2,
        "{modis_date1}": modis_date1,
        "{modis_date2}": modis_date2,
        "{n_cores}": n_cores,
        "{peak_memory}": peak_memory,
        "{run_name}": run_name,
        "{runtime}": f"{runtime:.2f}",
        "{spatial_param}": spatial_param,
        "{temporal_param}": temporal_param,
        "{tile_string}": tile_string,
        "{shapefile}": aoi_label,
    }

    # Read in the template
    with open(SUMMARY_TEMPLATE, "r") as template:
        lines = template.readlines()

    # Custom work for the files: only list outputs from this specific run
    lines_w_files = add_file_list(lines, output_directory, run_name, aoi,
                                  start_year, end_year,
                                  spatial_param, temporal_param)

    # Format template
    formatted_lines = []
    for line in lines_w_files:
        formatted_line = replace_values(parameters, line)
        formatted_lines.append(formatted_line)

    # Write to a README in the outputs directory
    sp = int(spatial_param)
    tp = int(temporal_param)
    aoi_part = f"_{aoi}" if aoi else ""
    fname = (f"{run_name.lower()}{aoi_part}_{start_year}-{end_year}"
             f"_s{sp:02d}_t{tp:02d}_readme.txt")
    fpath = output_directory.joinpath(fname)
    with open(fpath, "w") as summary:
        for line in formatted_lines:
            summary.write(line)
    logger.info(f"Run summary written to {fpath}.")
