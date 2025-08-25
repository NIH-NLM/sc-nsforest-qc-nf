# container/nsforest/context/src/nsforest_cli/violinplot.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import scanpy as sc
import matplotlib.pyplot as plt

from .utils import markers_from_nsforest_results
from .ensembl_lookup import ensg_to_symbol
from .dendro_subset import leaves_from_dendrogram


def violinplot_run(
    h5ad_in: Path,
    results_csv: Path,
    label_key: str,
    *,
    # existing
    markers_col: str = "NSForest_markers",
    cluster_col: str = "clusterName",
    clusters: Optional[Sequence[str]] = None,
    top_n: Optional[int] = None,
    log1p: bool = True,
    use_ensembl: bool = True,
    id_source: str = "var_names",
    dpi: int = 300,
    png_out: Optional[Path] = None,
    svg_out: Optional[Path] = None,
    # new leaf-driven subsetting
    leaf_range: Optional[str] = None,          # e.g., "0:10"
    leaf_indices: Optional[Sequence[int]] = None,  # e.g., 0 1 2 9
):
    adata = sc.read_h5ad(str(h5ad_in))

    # 1) Resolve clusters to keep
    if leaf_range or leaf_indices:
        # compute leaf-ordered labels and subset clusters by positions
        selected_labels = leaves_from_dendrogram(
            adata, label_key, leaf_range=leaf_range, leaf_indices=leaf_indices
        )
    elif clusters:
        selected_labels = list(clusters)
    else:
        selected_labels = None

    if selected_labels:
        adata = adata[adata.obs[label_key].isin(selected_labels)].copy()

    # 2) Markers → only for selected clusters (if any)
    markers_dict: Dict[str, List[str]] = markers_from_nsforest_results(
        results_csv, cluster_col=cluster_col, markers_col=markers_col
    )
    if not markers_dict:
        raise ValueError(f"No marker sets found in {results_csv}")

    if selected_labels:
        markers_dict = {k: v for k, v in markers_dict.items() if k in set(selected_labels)}

    # 3) Flatten marker gene list; de-dup; apply top_n per cluster
    gene_order: List[str] = []
    for cl, lst in markers_dict.items():
        take = lst if top_n is None else lst[:top_n]
        gene_order.extend(take)

    seen = set()
    gene_order = [g for g in gene_order if not (g in seen or seen.add(g))]
    present = [g for g in gene_order if g in adata.var_names]
    if not present:
        raise ValueError("None of the markers from results CSV are present in .var_names")

    adata = adata[:, present].copy()

    # 4) Symbols for display only
    display_names = list(adata.var_names)
    if use_ensembl:
        mapping = ensg_to_symbol(adata.var_names)
        display_names = [mapping.get(g, g) for g in adata.var_names]

    # 5) Log transform (as in NSForest tutorial)
    if log1p:
        sc.pp.log1p(adata)

    # 6) Ensure dendrogram is computed on the *subset* and show it in the violin
    sc.tl.dendrogram(adata, groupby=label_key)
    sc.pl.violin(
        adata,
        keys=adata.var_names.tolist(),
        groupby=label_key,
        log=False,
        use_raw=False,
        stripplot=False,
        dendrogram=True,
        rotation=90,
        show=False,
    )
    fig = plt.gcf()

    # 7) Swap x tick labels to symbols when counts match
    try:
        for ax in fig.axes:
            if len(ax.get_xticklabels()) == len(display_names):
                ax.set_xticklabels(display_names, rotation=90)
    except Exception:
        pass

    if png_out:
        fig.savefig(str(png_out), bbox_inches="tight", dpi=dpi)
    if svg_out:
        fig.savefig(str(svg_out), bbox_inches="tight", format="svg")

    return fig

