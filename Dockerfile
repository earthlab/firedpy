FROM continuumio/miniconda3
ENV PYTHONDONTWRITEBYTECODE=true

COPY . /home/firedpy

RUN rm -rf .git

WORKDIR /home/firedpy

RUN apt-get update && apt-get install -y --no-install-recommends screen htop curl unzip nano vim tree

# Download AWS CLI v2 and install it
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install

# Clean up the downloaded files and temporary packages
RUN rm -rf awscliv2.zip ./aws \
    && apt-get remove -y curl unzip \
    && apt-get clean

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

RUN conda update conda --yes \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict \
    && conda env create -f environment.yaml

RUN conda clean --all --yes --force-pkgs-dirs \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && conda list

# Set the PYTHONPATH environment variable
ENV PYTHONPATH="/home/firedpy"

RUN chmod +x /home/firedpy/entrypoint.sh
RUN echo "alias firedpy='python /home/firedpy/bin/firedpy.py'" >> /root/.bashrc
RUN echo "source activate firedpy" >> /root/.bashrc
# Use a custom entrypoint script to activate the conda environment
ENTRYPOINT ["/home/firedpy/entrypoint.sh"]
CMD ["/bin/bash"]
