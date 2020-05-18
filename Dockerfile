FROM continuumio/miniconda3:4.6.14

ENV PYTHONDONTWRITEBYTECODE=true

COPY environment.yaml .
COPY setup.py .

RUN conda update conda --yes \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict \
    && conda env create -f environment.yaml

RUN conda activate firedpy \
    && python setup.py install

RUN conda clean --all --yes --force-pkgs-dirs \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && conda list

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    awscli \
    htop
