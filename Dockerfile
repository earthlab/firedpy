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

# Download AWS CLI v2 and install it
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install

# Clean up the downloaded files and temporary packages
RUN rm -rf awscliv2.zip ./aws \
    && apt-get remove -y curl unzip \
    && apt-get clean

SHELL ["conda", "run", "-n", "firedpy", "/bin/bash", "-c"]

RUN python setup.py install 
