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
COPY context/requirements-freeze.txt .

# Install nsforest-cli CLI tool
RUN pip install ./

# Clone NSForest repo and expose it to Python
RUN git clone https://github.com/JCVenterInstitute/NSForest.git /opt/nsforest
ENV PYTHONPATH="/opt/nsforest/nsforest"

# Default: show CLI help
CMD ["nsforest-cli", "--help"]
