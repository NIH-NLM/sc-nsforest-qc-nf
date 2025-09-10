# nsforest_cli/dendrogramplot_run.py
from __future__ import annotations

from pathlib import Path
from typing import Optional, List

import matplotlib.pyplot as plt
import scanpy as sc
import nsforest as ns
import typer
from nsforest import as ns


def dendrogramplot_run(
        *,
        h5ad_in:       Path,
        label_key:     str,
        h5ad_out:      Path,
        png_out:       Optional[Path],
        svg_out:       Optional[Path],
        leaf_range:    Optional[str],
        leaf_indices:  Optional[List[int]],):
        
    """
    Plot dendrogramplot   
    """
    adata = sc.read_h5ad(str(h5ad_in))

    if label_key not in adata.obs:
        typer.echo(f"obs['{label_key}'] not found in AnnData.", err=True)
        raise typer.Exit(code=2)


    # Persist dendrogram structure in `uns` (no files written here)
    ns.pp.dendrogram(adata, cluster_header=label_key, save=True)

    # because we are using this as part of a workflow - we save the persisted dendrogram in a new h5ad
    adata.write_h5ad(str(h5ad_out))

    return None

