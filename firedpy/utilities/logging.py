"""firedpy Logging Utilities.

Run init_logger here with a target directory after creating a logging object
(`logger = logging.getLogger(__name__)`).
"""
import logging

from pathlib import Path


FORMAT = "%(levelname)s - %(asctime)s [%(filename)s:%(lineno)d] : %(message)s"
LOG_LEVELS = {
    "INFO": logging.INFO,
    "DEBUG": logging.DEBUG,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL
}


def init_logger(project_directory, log_level="INFO"):
    """Initiate logging to a file for a given firedpy project directory.

    Parameters
    ----------
    project_directory : str | pathlib.PosixPath
        The output directory for the firedpy run.
    log_level : str
        Level of logging detail to write to the log file.
    """
    # Setup firedpy logging
    filename = Path(project_directory).joinpath("logs/firedpy.log").absolute()
    filename = filename.expanduser()
    filename.parent.mkdir(exist_ok=True, parents=True)
    log_level = log_level.upper()
    logging.basicConfig(
        filename=filename,
        format=FORMAT,
        level=LOG_LEVELS[log_level],
        encoding="utf-8"
    )

    # Control third party logging
    logging.getLogger("earthaccess").setLevel(logging.WARNING)
