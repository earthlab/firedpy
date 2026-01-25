"""firedpy Command Line Interface (CLI)."""
import click
import time
import firedpy

from firedpy.help import HELP
from firedpy.run import run_firedpy, cleanup_intermediate_files
from firedpy.utilities.create_readme import make_read_me


CONTEXT_SETTINGS = {"max_content_width": 100, "terminal_width": 100}


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=firedpy.__version__)
@click.argument("project_directory", required=True)
@click.option("-t", "--tiles", required=True, help=HELP["tiles"])
@click.option("-y1", "--start_year", default=2000, help=HELP["year1"])
@click.option("-y2", "--end_year", default=2025, help=HELP["year2"])
@click.option("-d", "--daily", is_flag=True, help=HELP["daily"])
@click.option("-sp", "--spatial_param", default=8, help=HELP["sp"])
@click.option("-tp", "--temporal_param", default=3, help=HELP["tmp"])
@click.option("-shp", "--shape_file", default=None, help=HELP["shp"])
@click.option("-st", "--shape_type", default="gpkg", help=HELP["st"])
@click.option("-el", "--eco_region_level", default=1, help=HELP["eco_level"])
@click.option("-et", "--eco_region_type", default="na", help=HELP["eco_type"])
@click.option("-lt", "--land_cover_type", default=None, help=HELP["lc_type"])
@click.option("-f", "--full_csv", is_flag=True, help=HELP["full_csv"])
@click.option("-n", "--n_cores", default=0, help=HELP["n_cores"])
@click.option("-c", "--cleanup", is_flag=True, help=HELP["cleanup"])
def firedpy(project_directory, tiles, start_year, end_year, daily,
            spatial_param, temporal_param, shape_file, shape_type,
            eco_region_level, eco_region_type, land_cover_type, n_cores,
            full_csv, cleanup):
    """firedpy command line interface."""
    # Start the timer (seconds, not as helpful for prompted inputs)
    start = time.perf_counter()

    # Get the maximum resource use for this process (fix this)
    # initial_memory = peak_memory()

    # Covert tile list from string to Python list
    tiles = tiles.split()

    # Run with user parameters
    gdf, base_file_name = run_firedpy(
        out_dir=project_directory,
        tiles=tiles,
        start_year=start_year,
        end_year=end_year,
        daily=daily,
        spatial_param=spatial_param,
        temporal_param=temporal_param,
        shape_file=shape_file,
        shape_type=shape_type,
        eco_region_level=eco_region_level,
        eco_region_type=eco_region_type,
        land_cover_type=land_cover_type,
        full_csv=full_csv,
        n_cores=n_cores
    )

    # Done.
    end = time.perf_counter()
    seconds = end - start
    minutes = seconds / 60
    print(f"Job completed in {minutes:.2f} minutes")

    # peak_mem = (peak_memory() - initial_memory) / 1_024 ** 3
    # print(f"Peak memory usage: {peak_mem:.2f} GB")

    make_read_me(
        gdf=gdf,
        out_dir=project_directory,
        tile_name=None,
        file_base=base_file_name,
        input=1,
        daily=daily,
        spatial_param=spatial_param,
        temporal_param=temporal_param,
        shapefile=shape_file,
        shp_type=shape_type,
        job_time=seconds,
        n_cores=n_cores,
        # job_memory=peak_mem,
    )

    # Remove intermediate files if requested
    if cleanup:
        cleanup_intermediate_files(project_directory)
