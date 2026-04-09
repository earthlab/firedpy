"""firedpy Command Line Interface (CLI)."""
import click
import re

from click.core import ParameterSource
from pathlib import Path, PosixPath

import firedpy

from firedpy.help import ATTR_HELP, CLI_HELP
from firedpy.run import fired


CONTEXT_SETTINGS = {
    "max_content_width": 110,
    "terminal_width": 110,
    "show_default": True
}
USER_SOURCES = [
    ParameterSource.ENVIRONMENT,
    ParameterSource.COMMANDLINE
]


def clean_params(params):
    """Adjust parameters for proper data types and specific translations.

    Parameters
    ----------
    params : dict
        Parameter name-value pairs chosen from Firedpy CLI user input or by
        default.

    Returns
    -------
    dict : A copy of the original parameter dictionary with needed adjustments
        and translations.
    """
    clean_params = {}
    for key, value in params.items():
        # We have to use "None" instead of None for default value acceptance
        if value == "None":
            value = None

        # Expand paths
        if key in ["project_directory", "shape_file"]:
            if value:
                value = Path(value).absolute().expanduser()

        # Infer the project name if not provided
        if key == "project_name" and value is None:
            proj_dir = Path(params["project_directory"])
            value = proj_dir.absolute().expanduser().name

        clean_params[key] = value

    return clean_params


def confirm_parameters(params):
    """Print chosen parameters out and confirm if interactive.

    Parameters
    ----------
    params : dict
        Parameter name-value pairs chosen from Firedpy CLI user input or by
        default.

    Returns
    -------
    bool : Boolean indicating whether to continue with the firedpy run or to
        cancel (only for interactive mode).
    """
    # Collect user vs default parameter value pairs
    defaults = {}
    user_provided = {}
    for key, value in params.items():
        source = click.get_current_context().get_parameter_source(key)

        # Make these human readable with more context
        if source not in USER_SOURCES and key != "project_directory":
            defaults[key] = value
        else:
            user_provided[key] = value

    # Alert user to choices so they understand the defaults if used
    if user_provided:
        click.echo("Firedpy run submitted with the following user-supplied "
                   "values: ")
        for key, value in user_provided.items():
            helpful_print(key, value)

    if defaults:
        click.echo("Using the following default values: ")
        if key in ATTR_HELP:
            value = ATTR_HELP[key]
        for key, value in defaults.items():
            helpful_print(key, value)

    # If interactive, prompt user for confirmation
    proceed = True
    if params["interactive"]:
        proceed = click.confirm("Confirm parameters and proceed?")

    return proceed


def helpful_print(key, value):
    """Print key value pair so that they're easy to read & paste into a script.

    Parameters
    ----------
    key : str
        A string representing a firedpy parameter name.
    value : str
        A value associated with the firedpy parameter key.
    """
    # Give the user units for spatial and temporal values specifcally
    if key == "spatial_param":
        value = (
            f"{value}  # pixels (nominally ~{value * 463:,.0f} m but varies "
            "by location)"  # Let's calculate this value!
        )
        click.echo(f"    {key} = {value}")
    elif key == "temporal_param":
        value = f"{value}  # days"
        click.echo(f"    {key} = {value}")

    # Otherwise, add a comment if we've defined a hint for this value
    elif key in ATTR_HELP:
        if isinstance(value, (str, PosixPath)) and value != "None":
            click.echo(f"    {key} = '{value}'  # {ATTR_HELP[key][value]}")
        else:
            click.echo(f"    {key} = {value}  # {ATTR_HELP[key][value]}")

    # Otherwise, just print it out
    else:
        if isinstance(value, (str, PosixPath)) and value != "None":
            click.echo(f"    {key} = '{value}'")
        else:
            click.echo(f"    {key} = {value}")


def _prompts(ctx, _, interactive):
    """Callback click function that prompts the user for parameter values."""
    if interactive:
        click.echo(
            "Using interactive mode, please enter parameter values as "
            "prompted. Press enter to accept defaults if availabile (value in "
            "square brackets)\n"
        )
        for key, default in ctx.params.items():
            source = click.get_current_context().get_parameter_source(key)
            if source not in USER_SOURCES:
                # Skip eco_region_level when world ecoregions are selected
                if (key == "eco_region_level"
                        and ctx.params.get("eco_region_type") == "world"):
                    continue

                # Clean up the for-CLI help descriptions
                prompt = CLI_HELP[key].replace("\n", "")
                prompt = re.sub(" + ", " ", prompt).strip()

                # Prompt the user but adjust certain parameters for legibility
                value = click.prompt(prompt, default=default)
                if key == "project_directory" and value == ".":
                    value = Path(value).absolute().expanduser()

                click.echo(f"    {key}={value}")
                ctx.params[key] = value

        click.echo("")

    return interactive


@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.version_option(version=firedpy.__version__)
@click.option(
    "-i", "--interactive",
    is_flag=True,
    callback=_prompts,
    help=CLI_HELP["interactive"]
)
@click.option(
    "-p", "--project_directory",
    required=True,
    is_eager=True,
    default=None,
    help=CLI_HELP["interactive"]
)
@click.option(
    "-n", "--project_name",
    default="None",
    is_eager=True,
    help=CLI_HELP["project_name"]
)
@click.option(
    "-c", "--country",
    default="None",
    is_eager=True,
    help=CLI_HELP["country"]
)
@click.option(
    "-t", "--tiles",
    default="None",
    is_eager=True,
    help=CLI_HELP["tiles"]
)
@click.option(
    "-sf", "--shape_file",
    default="None",
    is_eager=True,
    help=CLI_HELP["shape_file"]
)
@click.option(
    "-y1", "--start_year",
    default=2000,
    is_eager=True,
    help=CLI_HELP["start_year"]
)
@click.option(
    "-y2", "--end_year",
    default=2025,
    is_eager=True,
    help=CLI_HELP["end_year"]
)
@click.option(
    "-sp", "--spatial_param",
    default=8,
    is_eager=True,
    help=CLI_HELP["spatial_param"]
)
@click.option(
    "-tp", "--temporal_param",
    default=3,
    is_eager=True,
    help=CLI_HELP["temporal_param"]
)
@click.option(
    "-d", "--daily",
    is_flag=True,
    is_eager=True,
    help=CLI_HELP["daily"]
)
@click.option(
    "-st", "--shape_type",
    default="gpkg",
    is_eager=True,
    help=CLI_HELP["shape_type"]
)
@click.option(
    "-et", "--eco_region_type",
    default="na",
    is_eager=True,
    help=CLI_HELP["eco_region_type"]
)
@click.option(
    "-el", "--eco_region_level",
    default=1,
    is_eager=True,
    help=CLI_HELP["eco_region_level"]
)
@click.option(
    "-lt", "--land_cover_type",
    default=1,
    is_eager=True,
    help=CLI_HELP["land_cover_type"]
)
@click.option(
    "-csv", "--csv_type",
    default="none",
    type=click.Choice(["full", "events", "none"], case_sensitive=False),
    is_eager=True,
    help=CLI_HELP["csv_type"]
)
@click.option(
    "-nc", "--n_cores",
    default=0,
    is_eager=True,
    help=CLI_HELP["n_cores"]
)
@click.option(
    "-cu", "--cleanup",
    is_flag=True,
    is_eager=True,
    help=CLI_HELP["cleanup"]
)
def cli(**params):
    """firedpy command line interface."""
    # Raise an error if parameters were supplied without the project directory
    projdir = params["project_directory"]
    if not projdir:
        raise ValueError("A value for the `project_directory` is required.")

    # Clean and confirm parameter values
    params = clean_params(params)
    confirmed = confirm_parameters(params)

    if confirmed:

        # The interactive option isn't in the main function
        del params["interactive"]

        # Run with user parameters
        _ = fired(**params)
