FROM continuumio/miniconda3:4.6.14

ENV PYTHONDONTWRITEBYTECODE=true

COPY . /home/firedpy

WORKDIR /home/firedpy

RUN conda update conda --yes \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict \
    && conda env create -f environment.yaml 

RUN conda clean --all --yes --force-pkgs-dirs \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && conda list

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    awscli \
    htop 
    
# The following line of code solved a problem that apparently is now not happening, and now this creates its own problem.
# If one is trying to do a docker build, and gets an error involving libffi.so.7, uncomment the following lines.
# RUN ln -s /opt/conda/envs/firedpy/lib/libffi.so.6 /opt/conda/envs/firedpy/lib/libffi.so.7 \
#  && pip install ipython

SHELL ["conda", "run", "-n", "firedpy", "/bin/bash", "-c"]

RUN python setup.py install 
