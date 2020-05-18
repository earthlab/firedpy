FROM continuumio/miniconda3:4.6.14

ENV PYTHONDONTWRITEBYTECODE=true

RUN conda update conda --yes \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict \
    && conda install --yes \
    python=3.7 \
    beautifulsoup4 \
    dask \
    descartes \
    geopandas \    
    gdal \
    lxml \
    matplotlib \
    netcdf4 \
    numpy \
    pandas \
    pyyaml \
    rasterio \
    scipy \ 
    shapely \
    toolz \ 
    tqdm \
    xarray=0.11.3 \
    && conda clean --all --yes --force-pkgs-dirs \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && conda list

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    awscli \
    htop 
