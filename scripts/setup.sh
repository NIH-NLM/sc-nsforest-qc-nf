!/usr/bin/bash
mamba install dask -y
mamba install anndata -y
mamba install plotly -y
mamba install zarr -y
mamba install h5py -y
mamba install gh -y
mamba install matplotlib -y
mamba uninstall numpy -y
mamba install "numpy<2"
mamba install pandas
mamba install scanpy
