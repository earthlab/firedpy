# Changelog

All notable changes to FIREDpy will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.0.1] - 2026
### Changed
- Regoranized code to allow for Python package installation
    - Moved all Python scripts to `firedpy` folder
        - This represents the most impactful set of changes
        - Scripts in `src` moved directly into `firedpy`.
            - `src/spatial` moved to `firedpy/utilities/spatial`.
        - Scripts in `bin` moved directly into `firedpy` and `firedpy.py` moved to `run.py` to avoid name conflicts with package.
        - Scripts in `utilities` moved to `firedpy/utilities`.
    - Moved code within scripts in general reorganization effort
        - Moved "HELP" string constants from `src/__init__.py` to dictionary to `help.py`.
    - Moved data in `ref` moved into `firedpy/data` and is handled as embedded package data.

- Began process of refactoring to comply with Python Enhancement Proposal (PEP) recommendations:
    - **Notable PEP recommendations include:**
        - *PEP-8: Style Guide for Python Code*
            - This is somewhat flexible, but is generally considered as a set of standards, is beneficial to readability, and is built into code linting checks embedded in continous integration and development (CI/CD).
            - https://peps.python.org/pep-0008/
        - *PEP-257: Docstring Conventions*
            - Provides a human readable description of function/class purpose, arguments, and outputs.
            - Will enable automated documentation with tools such as Sphinx (https://www.sphinx-doc.org/en/master/).
            - https://peps.python.org/pep-0257/
- Removed most type hinting and began inevitable debate over utility of practice in this context.
    - Reasons include:
        - Inconsistent use of hints in existing code base.
        - Lack of consistent, dedicated development team for proper implementation (code will be managed by Geography students)
        - Redundancy with docstring type declarations
        - Priority of simplicity and quick development time over benefits of IDE functionality and convenience.

### Added
- Packaging code added to `pyproject.toml` file, which is recommended practice as per PEP-621:
    - https://peps.python.org/pep-0621/
- Simple GDAL discovery routine held in `setup.py` script:
    - Called through `pyproject.toml`
    - Enables firedpy to use existing system installations of GDAL/OGR
    - Keeps Python environment lighter weight and reduces installation times
    - Provides an alternative to conda for users who prefer native Python package installations (e.g.,  Python Installs Packages or pip)
    - Is a necessary step for hosting firedpy on the Python Packaging Index (PyPI): https://pypi.org/


## [2.0.0] - 2023

### Added
- Temporal range: Nov 2000 - July 2024
- Parameters
    - CONUS
        - Temporal: 5 days
        - Spatial: 1 pixel
    - GLOBAL
        - Temporal: 5 days
        - Sptail: 1 pixel

### Changed
- Changed spatial and temporal aggregtation to 5 days, 1 pixels for all countries
- Bug fixes
- Optimized edge case algorithm
- New MODIS v6.1
- Extended time series

---

## [1.0.0] - 2022

### Added
- Initial FIREDpy release
- Temporal range: 2001 - 2021
- Parameters
    - CONUS
        - Temporal: 11 days
        - Spatial: 5 pixels
    - GLOBAL
        - Temporal: variable
        - Sptail: variable
- Validations: CONUS burned area comparison using MTBS and NIFC perimeters

---

## Template for Future Releases

### [X.Y.Z] - YYYY-MM-DD

#### Added
- New features and capabilities
- New data sources or processing methods
- New output formats or analysis tools

#### Changed
- Improvements to existing functionality
- Performance optimizations
- API updates (backward compatible)

#### Removed
- Features that have been removed

#### Fixed
- Bug fixes and error corrections
- Performance improvements
- Documentation updates

---

## Versioning Guidelines

### MAJOR (X.0.0)
- Core algorithm improvements (e.g., Köppen optimization)
- Breaking changes to API or data formats
- Major architectural changes

### MINOR (X.Y.0)
- Individual improvements (e.g., cFRP integration, speed integration)
- New features that are backward compatible
- New data sources or processing methods

### PATCH (X.Y.Z)
- Bug fixes and small performance updates
- Documentation improvements
- Code refactoring without functional changes

---

*This changelog is maintained as part of the FIREDpy development workflow. See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines and [development_workflow.md](development_workflow.md) for versioning policies.*
