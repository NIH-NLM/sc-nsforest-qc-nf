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

# Copy requirements
COPY context/requirements-freeze.txt .

# Install Python packages into base
#RUN pip install --upgrade pip \
# && pip install --no-cache-dir -r requirements-freeze.txt

# Clone NSForest repo (optional: remove if pip package covers all usage)
RUN git clone https://github.com/JCVenterInstitute/NSForest.git && \
    cd NSForest && \
    pip install -e .

# Set environment path to pick up CLI
ENV PATH="/opt/conda/bin:$PATH"

# Simple direct entrypoint
ENTRYPOINT [""]
