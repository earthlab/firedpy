"""firedpy Command Line Interface (CLI).

They want prompts. I want it to be optional. Here's an example of how to do
that:

https://gist.github.com/jzbruno/d167d256f6c426e2b6fdf2dbfc1ba3c9

"""
import click
import re
import time
import tracemalloc

from pathlib import Path

import firedpy

from firedpy.help import HELP
from firedpy.run import fired


CONTEXT_SETTINGS = {
    "max_content_width": 100,
    "terminal_width": 100,
    "show_default": True
}


def prompts(ctx, _, interactive):
    if interactive:
        click.echo(
            "Using interactive mode, please enter parameter values as "
            "prompted. Press enter to accept defaults if availabile (value in "
            "square brackets)"
        )
        for key, default in ctx.params.items():
            prompt = HELP[key].replace("\n", "")
            prompt = re.sub(" + ", " ", prompt).strip()
            value = click.prompt(prompt, default=default)
            if key == "project_directory" and value == ".":
                value = Path(value).absolute().expanduser()
            print(f"  {key}={value}")
            ctx.params[key] = value
    else:
        click.echo("Not using interactive mode.")
    return interactive


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=firedpy.__version__)
@click.option("-i", "--interactive", is_flag=True, callback=prompts,
              help=HELP["interactive"])
@click.option("-p", "--project_directory", default=None, required=True,
              is_eager=True, help=HELP["interactive"])
@click.option("-t", "--tiles", required=True, default="'h08v04 h09v04'",
              is_eager=True, help=HELP["tiles"])
@click.option("-y1", "--start_year", default=2000, is_eager=True,
              help=HELP["start_year"])
@click.option("-y2", "--end_year", default=2025, is_eager=True,
              help=HELP["end_year"])
@click.option("-d", "--daily", is_flag=True, is_eager=True,
              help=HELP["daily"])
@click.option("-sp", "--spatial_param", default=8, is_eager=True,
              help=HELP["spatial_param"])
@click.option("-tp", "--temporal_param", default=3, is_eager=True,
              help=HELP["temporal_param"])
@click.option("-st", "--shape_type", default="gpkg", is_eager=True,
              help=HELP["shape_type"])
@click.option("-el", "--eco_region_level", default=1, is_eager=True,
              help=HELP["eco_region_level"])
@click.option("-et", "--eco_region_type", default="na", is_eager=True,
              help=HELP["eco_region_type"])
@click.option("-lt", "--land_cover_type", default=1, is_eager=True,
              help=HELP["land_cover_type"])
@click.option("-f", "--full_csv", is_flag=True, is_eager=True,
              help=HELP["full_csv"])
@click.option("-n", "--n_cores", default=0, is_eager=True,
              help=HELP["n_cores"])
@click.option("-c", "--cleanup", is_flag=True, is_eager=True,
              help=HELP["cleanup"])
@click.pass_context
def cli(ctx, **params):
    """firedpy command line interface."""
    # Start the timer and memory tracer
    start = time.perf_counter()
    tracemalloc.start()

    # The interactive flag isn't present in the main function
    del params["interactive"]

    # Run with user parameters
    _ = fired(**params)

    # Get the maximum resource use for this process (fix this)
    _, peak_memory = tracemalloc.get_traced_memory()
    peak_memory = round(peak_memory / 1024 ** 3, 2)
    tracemalloc.stop()

    # Done.
    end = time.perf_counter()
    seconds = end - start
    runtime = seconds / 60
    print(f"Job completed in {runtime:.2f} minutes")
    print(f"Peak memory usage: {peak_memory:.2f} GB")
