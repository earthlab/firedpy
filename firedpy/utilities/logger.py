"""firedpy Logging Utilities.

Run init_logger here with a target directory after creating a logging object
(`logger = logging.getLogger(__name__)`).
"""
import logging
import os
import sys

from pathlib import Path


FORMAT = "%(levelname)s - %(asctime)s [%(filename)s:%(lineno)d] : %(message)s"
LOG_LEVELS = {
    "INFO": logging.INFO,
    "DEBUG": logging.DEBUG,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL
}


def init_logger(project_directory, run_name=None, aoi=None,
                start_year=None, end_year=None,
                spatial_param=None, temporal_param=None,
                log_level="INFO"):
    """Initiate logging to a file for a given firedpy project directory.

    Parameters
    ----------
    project_directory : str | pathlib.PosixPath
        The output directory for the firedpy run.
    run_name : str | None
        A name identifying this run, used in the log filename. If None,
        falls back to the project directory name.
    aoi : str | None
        Area-of-interest label (e.g. country name or shapefile
        stem), included in the log filename.
    start_year : int | None
        The first year of the run, included in the log filename.
    end_year : int | None
        The last year of the run, included in the log filename.
    spatial_param : int | None
        Spatial parameter for the run, included in the log filename.
    temporal_param : int | None
        Temporal parameter for the run, included in the log filename.
    log_level : str
        Level of logging detail to write to the log file.
    """
    # Build a descriptive log file name tied to this run
    if not run_name:
        run_name = Path(project_directory).name
    pid = os.getpid()
    if aoi and start_year and end_year and spatial_param and temporal_param:
        sp = int(spatial_param)
        tp = int(temporal_param)
        fname = (f"{run_name}_{aoi}_{start_year}-{end_year}"
                 f"_s{sp:02d}_t{tp:02d}_{pid}.log")
    else:
        fname = f"{run_name}_{pid}.log"
    fpath = Path(project_directory).joinpath("outputs", fname).absolute()
    fpath = fpath.expanduser()
    fpath.parent.mkdir(exist_ok=True, parents=True)

    # Make sure log level specific matches case in the level dictionary
    log_level = log_level.upper()

    # Setup firedpy logging
    logging.basicConfig(
        format=FORMAT,
        level=LOG_LEVELS[log_level],
        encoding="utf-8",
        handlers=[
            logging.FileHandler(fpath),
            logging.StreamHandler(sys.stdout)
        ]
    )

    logging.getLogger("earthaccess").setLevel(logging.WARNING)
