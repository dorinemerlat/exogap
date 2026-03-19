FROM mambaorg/micromamba:1.5.10

LABEL maintainer="Dorine MERLAT <dorine.merlat@etu.unistra.fr>"

USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        libgl1 \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        tcl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/GonzalezLab/MCHelper.git /opt/MCHelper && \
    chown -R $MAMBA_USER:$MAMBA_USER /opt/MCHelper

USER $MAMBA_USER

ENV MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba create -y -n MCHelper -f /opt/MCHelper/MCHelper.yml && \
    micromamba clean --all --yes

SHELL ["/usr/local/bin/_dockerfile_shell.sh"]

# Très important pour Singularity/runtime
ENV PATH=/opt/conda/envs/MCHelper/bin:/opt/conda/bin:$PATH
ENV CONDA_DEFAULT_ENV=MCHelper

CMD ["/bin/bash"]