FROM continuumio/miniconda3:22.11.1

ENV PYTHONDONTWRITEBYTECODE=true

COPY . /home/firedpy

WORKDIR /home/firedpy

RUN conda update conda --yes \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict \
    && conda env create -f environment.yaml \ 
    && echo "conda activate firedpy" >> ~/.bashrc

RUN conda clean --all --yes --force-pkgs-dirs \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && conda list

RUN apt-get update \
    && apt-get install -y htop curl unzip

# Download AWS CLI v2 and install it
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install

# Clean up the downloaded files and temporary packages
RUN rm -rf awscliv2.zip ./aws \
    && apt-get remove -y curl unzip \
    && apt-get clean

# The following line of code solved a problem that apparently is now not happening, and now this creates its own problem.
# If one is trying to do a docker build, and gets an error involving libffi.so.7, uncomment the following lines.
# RUN ln -s /opt/conda/envs/firedpy/lib/libffi.so.6 /opt/conda/envs/firedpy/lib/libffi.so.7 \
#  && pip install ipython

SHELL ["conda", "run", "-n", "firedpy", "/bin/bash", "-c"]

RUN python setup.py install  
