# violinplot.py

import scanpy as sc
import nsforest as ns
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from nsforest_cli.utils_load_convert_markers import load_and_convert_markers,convert_adata_varnames_with_symbol_map

def violinplot_run(
    h5ad_in,
    results_csv,
    symbol_map_csv=None,
    label_key="author_cell_type",
    svg_out=None,
    png_out=None,
    leaf_indices=None,
    leaf_range=None,
):
    # Load h5ad
    adata = sc.read(h5ad_in)

    # Load and convert markers
    markers_dict = load_and_convert_markers(results_csv, symbol_map_csv)

    # Convert var_names in adata if symbol_map provided
    if symbol_map_csv is not None:
        adata = convert_adata_varnames_with_symbol_map (adata, symbol_map_csv)
        adata.raw = None

    # Filter markers to those in adata.var_names (NOT raw.var_names)
    adata_genes = set(adata.var_names)
    filtered_markers = {
        k: [g for g in v if g in adata_genes]
        for k, v in markers_dict.items()
    }
    filtered_markers = {k: v for k, v in filtered_markers.items() if v}

    num_total = sum(len(v) for v in markers_dict.values())
    num_filtered = sum(len(v) for v in filtered_markers.values())
    print(f"[matrixplot] {num_filtered}/{num_total} markers found in adata.var_names")

    if not filtered_markers:
        raise ValueError("No marker genes found in AnnData. Check your symbol map or gene IDs.")

    # Persist dendrogram structure in `uns` (no files written here)
    fig = plt.figure()

    # Plot
    ax  = ns.pl.stackedviolin (
        adata,
        filtered_markers,
        cluster_header=label_key,
        show=False,
        save=True,
        output_folder = ".",
        outputfilename_suffix = label_key,
        use_raw=False,
    )

    # capture the figure that was actually drawn on
    if hasattr(ax, "get_figure"):
        fig = ax.get_figure()
    elif isinstance(ax, (list, tuple)) and ax and hasattr(ax[0], "get_figure"):
        fig = ax[0].get_figure()
    else:
        fig = plt.gcf()

    if png_out:
        fig.savefig(str(png_out), bbox_inches="tight")
    if svg_out:
        fig.savefig(str(svg_out), bbox_inches="tight", format="svg")

    plt.close(fig)

    return None
