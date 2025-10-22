# Changelog

All notable changes to FIREDpy will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

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
