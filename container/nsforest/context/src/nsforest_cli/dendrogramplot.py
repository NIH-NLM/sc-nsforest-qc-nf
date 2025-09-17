# nsforest_cli/dendrogramplot_run.py
from __future__ import annotations
import pandas as pd

from pathlib import Path
from typing import Optional, List

import matplotlib.pyplot as plt
import scanpy as sc
import nsforest as ns
import typer


def dendrogramplot_run(
        *,
        h5ad_in:        Path,
        label_key:      str,
        symbol_map_csv: Path,
        h5ad_out:       Path,
        leaf_range:     Optional[str],
        leaf_indices:   Optional[List[int]],):
        
    """
    Plot NSForest dendrogramplot   
    """

    # load h5ad
    adata = sc.read_h5ad(str(h5ad_in))

    # grab the base-prefix of the h5ad file
    # add the plot type here.
    base_prefix = h5ad_in.stem  
    suffix = "dendrogramplot"

    # Final output filename
    outputfilename_suffix = f"{base_prefix}-{suffix}"

    # Load and check symbol map CSV
    symbol_map_df = pd.read_csv(symbol_map_csv)

    required_columns = {'ensg', 'symbol'}
    if not required_columns.issubset(symbol_map_df.columns):
        typer.echo(f"[ERROR] Symbol map must have columns: {required_columns}", err=True)
        raise typer.Exit(code=2)

    symbol_map = dict(zip(symbol_map_df['ensg'], symbol_map_df['symbol']))

    # Backup original var_names
    adata.var["orig_names"] = adata.var_names

    # Safe mapping: replace Ensembl with gene symbols if possible
    adata.var_names = [
        symbol_map.get(gene_id, gene_id) for gene_id in adata.var_names
    ]

    # Warn if nothing mapped
    num_mapped = sum(1 for orig, new in zip(adata.var["orig_names"], adata.var_names) if orig != new)

    if num_mapped == 0:
        typer.echo("[WARNING] No Ensembl IDs were mapped to gene symbols. Check your symbol map CSV.", err=True)
    
    if label_key not in adata.obs:
        typer.echo(f"obs['{label_key}'] not found in AnnData.", err=True)
        raise typer.Exit(code=2)


    # Persist dendrogram structure in `uns` (no files written here)
    ns.pp.dendrogram(
        adata,
        cluster_header=label_key,
        save=True,
        output_folder = ".",
        outputfilename_suffix= outputfilename_suffix)

    # because we are using this as part of a workflow - we save the persisted dendrogram in a new h5ad
    adata.write_h5ad(str(h5ad_out))

    return None

