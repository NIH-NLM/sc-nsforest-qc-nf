# container/nsforest/context/src/nsforest_cli/dendrogram.py
from __future__ import annotations

from pathlib import Path
from typing import Optional
import scanpy as sc
import matplotlib.pyplot as plt
from .dendro_subset import leaves_from_dendrogram


def dendrogramplot_run(
    h5ad_in: Path,
    results_csv: Path,   # kept for uniform CLI; not used here
    label_key: str,
    *,
    png_out: Optional[Path] = None,
    svg_out: Optional[Path] = None,
    dpi: int = 300,
    # new leaf subsetting
    leaf_range: Optional[str] = None,
    leaf_indices: Optional[list[int]] = None,
):
    adata = sc.read_h5ad(str(h5ad_in))

    # optional subsetting by leaf positions
    if leaf_range or leaf_indices:
        selected_labels = leaves_from_dendrogram(
            adata, label_key, leaf_range=leaf_range, leaf_indices=leaf_indices
        )
        adata = adata[adata.obs[label_key].isin(selected_labels)].copy()

    sc.tl.dendrogram(adata, groupby=label_key)
    sc.pl.dendrogram(adata, groupby=label_key, show=False)
    fig = plt.gcf()

    if png_out:
        fig.savefig(str(png_out), bbox_inches="tight", dpi=dpi)
    if svg_out:
        fig.savefig(str(svg_out), bbox_inches="tight", format="svg")

    return fig

