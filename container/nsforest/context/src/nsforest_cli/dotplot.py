# container/nsforest/context/src/nsforest_cli/dotplot.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence
import scanpy as sc
import matplotlib.pyplot as plt

from .utils import markers_from_nsforest_results
from .ensembl_lookup import ensg_to_symbol
from .dendro_subset import leaves_from_dendrogram


def dotplot_run(
    h5ad_in: Path,
    results_csv: Path,
    label_key: str,
    *,
    markers_col: str = "NSForest_markers",
    cluster_col: str = "clusterName",
    clusters: Optional[Sequence[str]] = None,
    top_n: Optional[int] = None,
    log1p: bool = True,
    use_ensembl: bool = True,
    dpi: int = 300,
    png_out: Optional[Path] = None,
    svg_out: Optional[Path] = None,
    # new leaf subsetting
    leaf_range: Optional[str] = None,
    leaf_indices: Optional[Sequence[int]] = None,
):
    adata = sc.read_h5ad(str(h5ad_in))

    # 1) cluster selection
    if leaf_range or leaf_indices:
        selected_labels = leaves_from_dendrogram(
            adata, label_key, leaf_range=leaf_range, leaf_indices=leaf_indices
        )
    elif clusters:
        selected_labels = list(clusters)
    else:
        selected_labels = None

    if selected_labels:
        adata = adata[adata.obs[label_key].isin(selected_labels)].copy()

    # 2) markers
    markers_dict: Dict[str, List[str]] = markers_from_nsforest_results(
        results_csv, cluster_col=cluster_col, markers_col=markers_col
    )
    if not markers_dict:
        raise ValueError(f"No marker sets found in {results_csv}")

    if selected_labels:
        markers_dict = {k: v for k, v in markers_dict.items() if k in set(selected_labels)}

    # 3) gene list
    gene_order: List[str] = []
    for _, lst in markers_dict.items():
        take = lst if top_n is None else lst[:top_n]
        gene_order.extend(take)
    seen = set()
    gene_order = [g for g in gene_order if not (g in seen or seen.add(g))]
    present = [g for g in gene_order if g in adata.var_names]
    if not present:
        raise ValueError("None of the markers from results CSV are present in .var_names")

    adata = adata[:, present].copy()

    # 4) log
    if log1p:
        sc.pp.log1p(adata)

    # 5) symbols for y ticklabels
    display_names = None
    if use_ensembl:
        mapping = ensg_to_symbol(adata.var_names)
        display_names = [mapping.get(g, g) for g in adata.var_names]

    dp = sc.pl.dotplot(
        adata,
        var_names=adata.var_names.tolist(),
        groupby=label_key,
        standard_scale="var",
        show=False,
        return_fig=True,
    )
    fig = dp.fig

    if use_ensembl and display_names:
        try:
            for ax in fig.axes:
                labs = [t.get_text() for t in ax.get_yticklabels()]
                if set(labs) == set(adata.var_names):
                    ax.set_yticklabels(display_names)
        except Exception:
            pass

    if png_out:
        fig.savefig(str(png_out), bbox_inches="tight", dpi=dpi)
    if svg_out:
        fig.savefig(str(svg_out), bbox_inches="tight", format="svg")

    return fig

