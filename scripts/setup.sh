!/mnt/libraries/envs/cloudos/bin/bash
mamba install emacs -y
mamba install dask -y
mamba install anndata -y
mamba install plotly -y
mamba install zarr -y
mamba install h5py -y
mamba install gh -y
mamba install matplotlib -y
mamba uninstall numpy -y
mamba install "numpy<2" -y
mamba install pandas -y
mamba install scanpy -y
