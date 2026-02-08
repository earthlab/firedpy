"""firedpy Command Line Interface (CLI).

They want prompts. I want it to be optional. Here's an example of how to do
that:

https://gist.github.com/jzbruno/d167d256f6c426e2b6fdf2dbfc1ba3c9

"""
import click
import os
import re
import time
import tracemalloc

from click.core import ParameterSource
from pathlib import Path

import firedpy

from firedpy.help import ATTR_HELP, CLI_HELP
from firedpy.run import fired


CONTEXT_SETTINGS = {
    "max_content_width": 100,
    "terminal_width": 100,
    "show_default": True
}
USER_SOURCES = [
    ParameterSource.ENVIRONMENT,
    ParameterSource.COMMANDLINE
]


def print_parameters(params):
    """Print the user and default parameters out from CLI context object."""
    # Collect user vs default parameter value pairs
    defaults = {}
    user_provided = {}
    for key, value in params.items():
        source = click.get_current_context().get_parameter_source(key)

        # Make these human readable with more context
        if key in ATTR_HELP:
            value = ATTR_HELP[key][value]
        if key == "temporal_param":
            value = f"{value} days"
        if key == "spatial_param":
            value = (
                f"{value} pixels. Nominally ~{value * 463:,.0f} m but varies "
                "by location)"  # Let's calculate this value!
            )
        if source not in USER_SOURCES:
            defaults[key] = value
        else:
            user_provided[key] = value

    # Alert user to choices so they understand the defaults if used
    if user_provided:
        click.echo("Firedpy run submitted with the following user-supplied "
                   "values: ")
        for key, value in user_provided.items():
            click.echo(f"    {key}: {value}")
    if defaults:
        click.echo("Using the following default values: ")
        if key in ATTR_HELP:
            value = ATTR_HELP[key]
        for key, value in defaults.items():
            click.echo(f"    {key}: {value}")

    # Give them a few seconds to decide to cancel the process
    wait = 3
    click.echo(f"Running firedpy in {wait} seconds...")
    time.sleep(wait)


def _prompts(ctx, _, interactive):
    """Callback click function that prompts the user for parameter values."""
    if interactive:
        click.echo(
            "Using interactive mode, please enter parameter values as "
            "prompted. Press enter to accept defaults if availabile (value in "
            "square brackets)"
        )
        for key, default in ctx.params.items():
            source = click.get_current_context().get_parameter_source(key)
            if source not in USER_SOURCES:
                prompt = CLI_HELP[key].replace("\n", "")
                prompt = re.sub(" + ", " ", prompt).strip()
                value = click.prompt(prompt, default=default)
                if key == "project_directory" and value == ".":
                    value = Path(value).absolute().expanduser()
                if key == "n_cores" and value == 0:
                    value = os.cpu_count()
                click.echo(f"  {key}={value}")
                ctx.params[key] = value
    return interactive


@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.version_option(version=firedpy.__version__)
@click.option("-i", "--interactive", is_flag=True, callback=_prompts,
              help=CLI_HELP["interactive"])
@click.option("-p", "--project_directory", default=None, required=True,
              is_eager=True, help=CLI_HELP["interactive"])
@click.option("-t", "--tiles", required=True, default="'h08v04 h09v04'",
              is_eager=True, help=CLI_HELP["tiles"])
@click.option("-y1", "--start_year", default=2000, is_eager=True,
              help=CLI_HELP["start_year"])
@click.option("-y2", "--end_year", default=2025, is_eager=True,
              help=CLI_HELP["end_year"])
@click.option("-d", "--daily", is_flag=True, is_eager=True,
              help=CLI_HELP["daily"])
@click.option("-sp", "--spatial_param", default=8, is_eager=True,
              help=CLI_HELP["spatial_param"])
@click.option("-tp", "--temporal_param", default=3, is_eager=True,
              help=CLI_HELP["temporal_param"])
@click.option("-st", "--shape_type", default="gpkg", is_eager=True,
              help=CLI_HELP["shape_type"])
@click.option("-el", "--eco_region_level", default=1, is_eager=True,
              help=CLI_HELP["eco_region_level"])
@click.option("-et", "--eco_region_type", default="na", is_eager=True,
              help=CLI_HELP["eco_region_type"])
@click.option("-lt", "--land_cover_type", default=1, is_eager=True,
              help=CLI_HELP["land_cover_type"])
@click.option("-f", "--full_csv", is_flag=True, is_eager=True,
              help=CLI_HELP["full_csv"])
@click.option("-n", "--n_cores", default=0, is_eager=True,
              help=CLI_HELP["n_cores"])
@click.option("-c", "--cleanup", is_flag=True, is_eager=True,
              help=CLI_HELP["cleanup"])
def cli(**params):
    """firedpy command line interface."""
    # Start the timer and memory tracer
    start = time.perf_counter()
    tracemalloc.start()

    # Raise an error if parameters were supplied without the project directory
    projdir = params["project_directory"]
    if not projdir:
        raise ValueError("A value for the `project_directory` is required.")

    # Make sure the project directory looks good if interactive wasn't used
    if projdir == "." or "~" in str(projdir):
        params["project_directory"] = Path(projdir).absolute().expanduser()

    # Make sure the cpu core count looks good if interactive wasn't used
    n_cores = params["n_cores"]
    if n_cores == 0:
        params["n_cores"] = os.cpu_count()

    # Print out the parameter values for non-user supplied defaults
    print_parameters(params)

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
    click.echo(f"Job completed in {runtime:.2f} minutes")
    click.echo(f"Peak memory usage: {peak_memory:.2f} GB")
