# container/nsforest/context/src/nsforest_cli/dendrogramplot.py
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import nsforest as ns
import scanpy as sc
import matplotlib.pyplot as plt


def _slice_leaves(
    leaves: list[str],
    *,
    leaf_range: Optional[str] = None,
    leaf_indices: Optional[Sequence[int]] = None,
) -> list[str]:
    """Subset dendrogram leaves by indices or by slice string 'start:end' (end-exclusive).
    If both provided, `leaf_indices` wins.
    """
    if leaf_indices:
        n = len(leaves)
        idx = [i for i in leaf_indices if -n <= i < n]
        return [leaves[i] for i in idx]
    if leaf_range:
        try:
            start_s, end_s = leaf_range.split(":", 1)
            start = int(start_s) if start_s.strip() != "" else None
            end = int(end_s) if end_s.strip() != "" else None
            return leaves[slice(start, end)]
        except Exception:
            return leaves
    return leaves


def dendrogramplot_run(
    h5ad_in: Path,
    label_key: str,
    png_out: str,
    svg_out: str,
    leaf_range: Optional[str] = None,
    leaf_indices: Optional[Sequence[int]] = None,
) -> None:
    
    # Read and validate
    adata = sc.read_h5ad(str(h5ad_in))
    if label_key not in adata.obs:
        raise ValueError(f"obs['{label_key}'] not found in {h5ad_in}")
    adata.obs[label_key] = adata.obs[label_key].astype("category")

    # Ensure we have a dendrogram on the full set to derive the leaf order
    dendro_key = f"dendrogram_{label_key}"
    if dendro_key not in adata.uns or "categories_ordered" not in adata.uns.get(dendro_key, {}):
        adata = ns.pp.dendrogram(adata,
                         cluster_header=label_key,
                         save = "svg",
                         output_folder = ".",
                         outputfilename_suffix = label_key)
        
    leaves = list(adata.uns[dendro_key]["categories_ordered"])

    # Apply optional leaf slicing
#    sliced = _slice_leaves(leaves, leaf_range=leaf_range, leaf_indices=leaf_indices)
#    selected = sliced

    # Need at least two categories for a dendrogram plot
#    if len(selected) < 2:
#        raise ValueError("Need at least two clusters to plot a dendrogram after subsetting.")

#   adata.obs[label_key] = adata.obs[label_key].cat.set_categories(selected, ordered=True)

    # Recompute dendrogram on the subset for a correct tree and plot
    sc.tl.dendrogram(adata, groupby=label_key)
    sc.pl.dendrogram(adata, groupby=label_key, show=False)

    # Save outputs
    fig = plt.gcf()
    if png_out:
        Path(png_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(png_out), bbox_inches="tight")
    if svg_out:
        Path(svg_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(svg_out), bbox_inches="tight", format="svg")

    plt.close(fig)
    return adata

