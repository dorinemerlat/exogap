FROM continuumio/miniconda3

LABEL maintainer="Dorine MERLAT <dorine.merlat@etu.unistra.fr>"

WORKDIR /opt

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        libgl1 \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/GonzalezLab/MCHelper.git /opt/MCHelper

RUN conda env create -f /opt/MCHelper/MCHelper.yml && \
    conda clean -afy

SHELL ["conda", "run", "-n", "MCHelper", "/bin/bash", "-c"]

RUN python -c "import cv2; print('cv2 OK')"

RUN printf '#!/bin/sh\nexec conda run --no-capture-output -n MCHelper python /opt/MCHelper/MCHelper.py \"$@\"\n' \
    > /usr/local/bin/MCHelper.py && \
    chmod +x /usr/local/bin/MCHelper.py

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "MCHelper", "python", "/opt/MCHelper/MCHelper.py"]