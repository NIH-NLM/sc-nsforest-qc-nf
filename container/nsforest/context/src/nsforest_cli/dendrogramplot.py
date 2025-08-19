# nsforest_cli/dendrogramplot.py
from pathlib import Path
from typing import Optional

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

def dendrogramplot_run(
    h5ad_in: Path,
    results_csv: Path,              # kept for interface consistency; not used currently
    label_key: str,
    png_out: Optional[Path] = None,
    svg_out: Optional[Path] = None,
):
    """
    Compute dendrogram (sc.tl.dendrogram) and plot it. Optionally save PNG/SVG.
    Returns the AnnData (so caller can write .h5ad).
    """
    adata = sc.read_h5ad(str(h5ad_in))

    # Ensure groupby column is categorical
    if not pd.api.types.is_categorical_dtype(adata.obs[label_key]):
        adata.obs[label_key] = adata.obs[label_key].astype("category")

    # Compute dendrogram (stores linkage & info in adata.uns[f"dendrogram_{label_key}"])
    sc.tl.dendrogram(adata, groupby=label_key)

    # Plot without trying to use return_fig (not available on all versions)
    ax = sc.pl.dendrogram(adata, groupby=label_key, show=False)
    # Resolve figure handle robustly
    if ax is None:
        fig = plt.gcf()
    else:
        fig = getattr(ax, "figure", plt.gcf())

    # Save, if requested
    if png_out:
        fig.savefig(str(png_out), bbox_inches="tight", dpi=300)
    if svg_out:
        fig.savefig(str(svg_out), bbox_inches="tight", format="svg")

    plt.close(fig)
    return adata

