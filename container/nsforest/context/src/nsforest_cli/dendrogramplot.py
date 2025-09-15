# nsforest_cli/dendrogramplot_run.py
from __future__ import annotations

from pathlib import Path
from typing import Optional, List

import matplotlib.pyplot as plt
import scanpy as sc
import nsforest as ns
import typer


def dendrogramplot_run(
        *,
        h5ad_in:       Path,
        label_key:     str,
        h5ad_out:      Path,
        leaf_range:    Optional[str],
        leaf_indices:  Optional[List[int]],):
        
    """
    Plot NSForest dendrogramplot   
    """

    # load h5ad
    adata = sc.read_h5ad(str(h5ad_in))

    # grab the base-prefix of the h5ad file
    # add the plot type here.
    base_prefix = h5ad_path.stem  
    suffix = "dendrogramplot"

    # Final output filename
    outputfilename_suffix = f"{base_prefix}-{suffix}"
    
    # Load the mapping file: assumes two columns: 'ensembl_id', 'gene_symbol'
    symbol_map_df = pd.read_csv("gencode-release-49-ensg-gene-symbol.csv")

    # Create a dict: ensembl_id → gene_symbol
    symbol_map = dict(zip(symbol_map_df['ensembl_id'], symbol_map_df['gene_symbol']))

    # Store original Ensembl IDs
    adata.var["orig_names"] = adata.var_names

    # Map var_names (which are currently Ensembl IDs) to gene symbols
    adata.var_names = adata.var_names.map(symbol_map).fillna(adata.var_names)
    
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

