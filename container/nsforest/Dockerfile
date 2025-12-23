# File: container/nsforest/Dockerfile

FROM mambaorg/micromamba:1.5.6

LABEL maintainer="nih-nlm"

USER root:root

RUN apt-get update && \
    apt-get install -y git procps && \
    apt-get install -y --no-install-recommends build-essential gcc g++ gfortran && \
    apt-get clean

WORKDIR /app

COPY context/nsforest.yml /app/env.yaml
COPY context/requirements-freeze.txt /app/requirements.txt
RUN chown -R mambauser:mambauser /app

COPY context/src ./src
RUN chown -R mambauser:mambauser ./src

RUN git clone https://github.com/JCVenterInstitute/NSForest.git /opt/nsforest && \
    chown -R mambauser:mambauser /opt/nsforest

USER mambauser:mambauser

ENV MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:$PATH \
    DEBIAN_FRONTEND=noninteractive

RUN micromamba install -y -n base -f /app/env.yaml && \
    micromamba clean --all --yes

COPY context/setup.py .

RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel
RUN python -m pip install --no-cache-dir -r /app/requirements.txt
RUN python -m pip install --no-cache-dir "plotly==5.22.0"
RUN python -m pip install --no-cache-dir /app/

ENV PYTHONPATH="/opt/nsforest/nsforest"

CMD ["nsforest-cli", "--help"]
