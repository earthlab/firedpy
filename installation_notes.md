### Implementing Python Installation Notes
- We're adding a pyproject.toml file because it's recommended now, but it looks like we'll still need to rely on a setup.py file for any dynamic installation steps.
- Building in Python 3.14 to catch up, but this will require a few dependency and code changes here and there.
- Moving all Python material and data into the "firedpy" folder.
- Would like to include a more standard installation guide in addition to conda. This should include the installation of any system level dependencies 
- This will need some more dependency version checks, especially since 3.14.
- Removed the minor release constraint on our dependencies, assuming we plan to continue maintaining it. 
- The firedpy version is a dynamically updated from its __init__.py
- Requirements are handled by setup.py, everything else by pyproject.toml
    - This allows us to match the system-level GDAL dependency so we're not fully dependent on conda.

