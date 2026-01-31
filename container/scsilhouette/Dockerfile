# File: container/scsilhouette/Dockerfile
FROM mambaorg/micromamba:latest

LABEL maintainer="nih-nlm"

USER root:root

RUN apt-get update && \
    apt-get install -y git procps build-essential gcc g++ gfortran python3-pip && \
    apt-get clean

WORKDIR /app

RUN git clone https://github.com/NIH-NLM/scsilhouette.git && \
    micromamba install -y -n base -f scsilhouette/environment.yml && \
    micromamba clean --all --yes && \
    cd scsilhouette && \
    micromamba run -n base pip install -e . && \
    chown -R mambauser:mambauser /app/scsilhouette

USER mambauser:mambauser

ENV PATH=/opt/conda/bin:$PATH \
    PYTHONPATH="/app/scsilhouette/src"

ENTRYPOINT ["scsilhouette"]
CMD ["--help"]
