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

RUN conda update conda --yes \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict \
    && conda env create -f environment.yml

RUN conda clean --all --yes --force-pkgs-dirs \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && conda list

# Set the PYTHONPATH environment variable
ENV PYTHONPATH="/home/firedpy"

RUN adduser --disabled-password --gecos "VICE_User" --uid 1000 user  && \
    usermod -aG sudo user && \
    echo "$LOCAL_USER ALL=NOPASSWD: $PRIV_CMDS" >> /etc/sudoers

RUN apt-get update && \
    apt-get install -y curl grep sed dpkg wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    gettext-base git mercurial subversion \
    tmux && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# install ttyd
RUN curl -L "https://github.com/tsl0922/ttyd/releases/download/1.6.3/ttyd.x86_64" > ttyd && \
    chmod a+x ttyd && \
    mv ttyd /usr/local/bin/ttyd

RUN apt-get update && \
    curl -L "https://github.com/krallin/tini/releases/download/v0.19.0/tini_0.19.0-amd64.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# set shell as bash and terminal as linux
ENV SHELL=bash
ENV TERM=xterm

# open port 7681 for ttyd
EXPOSE 7681

# changes tmux layout while running
COPY entry.sh /bin
RUN echo 'set-option -g status off' >> ~/.tmux.conf

RUN chmod +x /bin/entry.sh
RUN echo "alias firedpy='python /home/firedpy/bin/firedpy.py'" >> /root/.bashrc
RUN echo "source activate firedpy" >> /root/.bashrc
# Use a custom entrypoint script to activate the conda environment
ENTRYPOINT ["/bin/entry.sh"]
CMD ["ttyd", "bash"]
