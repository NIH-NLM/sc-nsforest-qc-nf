# container/nsforest/context/src/nsforest_cli/dotplot.py
from __future__ import annotations
import ast
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import nsforest as ns
import numpy as np
import scanpy as sc
from scipy import sparse
from typing import Dict, List, Optional, Sequence

def dotplot_run(
        *,
        h5ad_in: Path,
        results_csv: Path,
        label_key: str,
        png_out: Optional[Path] = None,
        svg_out: Optional[Path] = None,
        leaf_range: Optional[str] = None,
        leaf_indices: Optional[List[int]] = None,
) -> None:
    """
    Render a dotplot replicating the NSForest tutorial behavior.
    """

    # Load AnnData and ensure groupby column is categorical
    adata = sc.read_h5ad(str(h5ad_in))

    # Read NSForest results and align to cluster order
    df = pd.read_csv(results_csv)

    cluster_header = label_key
    dendrogram = list(adata.uns["dendrogram_" + cluster_header]["categories_ordered"])

    # Prepare to_plot DataFrame
    to_plot = df.copy()
    if "clusterName" in to_plot.columns:
        to_plot["clusterName"] = to_plot["clusterName"].astype("category")
        to_plot["clusterName"] = to_plot["clusterName"].cat.set_categories(dendrogram)
        to_plot = to_plot.sort_values("clusterName")
    if "NSForest_markers" in to_plot.columns:
        to_plot = to_plot.rename(columns={"NSForest_markers": "markers"})

    # Prepare markers_dict
    markers_dict = dict(zip(to_plot["clusterName"], to_plot["markers"]))

    ad_for_plot = adata
    save = True  

    # create & set the current figure
    fig = plt.figure()
    ax = ns.pl.dotplot(
        ad_for_plot,
        markers_dict,
        label_key,
        dendrogram=dendrogram,
        save=save,
        output_folder=".",
        outputfilename_suffix=".",
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

