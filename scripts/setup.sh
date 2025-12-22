#!/usr/bin/env bash

set -euo pipefail

ENV_NAME="nsforest-test"

echo "Creating environment: $ENV_NAME"

mamba create -n "$ENV_NAME" -c conda-forge \
  python=3.10 \
  scanpy \
  dask=2024.1.0 \
  numpy=1.26.4 \
  anndata \
  pandas \
  zarr \
  h5py \
  matplotlib -y

echo
echo "To activate the environment, run:"
echo "conda activate $ENV_NAME"
