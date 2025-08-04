# File: container/nsforest/Dockerfile

FROM mambaorg/micromamba:1.5.6

LABEL maintainer="nih-nlm"

ENV MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:$PATH \
    DEBIAN_FRONTEND=noninteractive

WORKDIR /app

# Copy environment file and install into base environment
COPY --chown=$MAMBA_USER:$MAMBA_USER context/nsforest.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Enable conda activation
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN echo "micromamba activate base" >> ~/.bashrc

# Temporarily become root to install OS packages
USER root
RUN apt-get update && \
    apt-get install -y git procps && \
    apt-get clean -y

# Copy source + setup
COPY context/ /app

# Install everything
RUN pip install ./

# Copy Python source
# This copies src to the current working directory /app/
# /app/
# |--- src/
# |     |__ nsforest_cli/
# |
# |--- requirements-freeze.txt
#
COPY context/src ./src

# PYTHONPATH is an environment variable that tells Python where to look
#  for importable packages and modules in addition to default system paths
#
# setting PYTHONPATH=/app/src is necessary for python -m nsforest_cli.main
#
ENV PYTHONPATH=/app/src


# Clone NSForest repo (optional: remove if pip package covers all usage)
RUN git clone https://github.com/JCVenterInstitute/NSForest.git /opt/nsforest && \
    pip install -e .

# Set environment path to pick up CLI
ENV PATH="/opt/conda/bin:$PATH"

# Simple direct entrypoint
ENTRYPOINT [""]

#
# We then can run from the command line in bash
#
# docker run --rm nsforest-cli prep-medians