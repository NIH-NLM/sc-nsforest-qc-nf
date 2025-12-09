# File: container/nsforest/Dockerfile

FROM mambaorg/micromamba:1.5.6

LABEL maintainer="nih-nlm"

ENV MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:$PATH \
    DEBIAN_FRONTEND=noninteractive

WORKDIR /app

# Install Conda environment
COPY context/nsforest.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

USER root
RUN apt-get update && \
    apt-get install -y git procps && \
    apt-get clean

# Copy CLI code
COPY context/src ./src
COPY context/setup.py .

# system toolchain for building wheels
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gcc g++ gfortran \
 && rm -rf /var/lib/apt/lists/*

# modern build frontends help many packages
RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel

# now install your requirements
COPY context/requirements-freeze.txt /app/requirements.txt
RUN python -m pip install --no-cache-dir -r /app/requirements.txt

# and force the Plotly you want (ensures from_matplotlib exists)
RUN python -m pip install --no-cache-dir "plotly==5.22.0"

# Install nsforest-cli CLI tool
RUN pip install /app/

# Clone NSForest repo and expose it to Python
RUN git clone https://github.com/adeslatt/NSForest.git /opt/nsforest && \
    cd /opt/nsforest && \
    git checkout 6cbe49b6c37f5921db87542d8ee53728370a9dbb

ENV PYTHONPATH="/opt/nsforest/nsforest"


# Default: show CLI help
CMD ["nsforest-cli", "--help"]
