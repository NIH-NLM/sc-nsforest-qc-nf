# container/nsforest/context/src/nsforest_cli/violinplot.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import scanpy as sc
import matplotlib.pyplot as plt


def violinplot_run(
    h5ad_in: Path,
    results_csv: Path,
    *,
    label_key: str,
    png_out: Optional[Path] = None,
    svg_out: Optional[Path] = None,
) -> None:
    """
    Render a stacked violin plot replicating the NSForest tutorial behavior.

    Behavior:
      - Reads an AnnData .h5ad and the NSForest results .csv.
      - Ensures a dendrogram exists for `label_key` and respects its category order.
      - Builds a markers dictionary from the results CSV.
      - Uses `nsforest.pl.stackedviolin(..., save=False)` to obtain a Matplotlib figure.
      - Saves exactly the filenames provided by `--png-out` / `--svg-out`.
      - Closes the figure and returns None.
    """
    import pandas as pd
    import ast
    import numpy as np
    from scipy import sparse
    import nsforest as ns

    # Load AnnData and ensure groupby column is categorical
    adata = sc.read_h5ad(str(h5ad_in))
    if label_key not in adata.obs:
        raise ValueError(f"obs['{label_key}'] not found in {h5ad_in}")
    adata.obs[label_key] = adata.obs[label_key].astype("category")

    # Ensure dendrogram exists; derive order
    dendro_key = f"dendrogram_{label_key}"
    if dendro_key not in adata.uns or "categories_ordered" not in adata.uns.get(dendro_key, {}):
        sc.tl.dendrogram(adata, groupby=label_key)
    dendrogram: List[str] = list(adata.uns[dendro_key]["categories_ordered"])

    # Optional subset of clusters while preserving dendrogram order
#    if clusters:
#        keep = set(clusters)
#        dendrogram = [c for c in dendrogram if c in keep]
#        adata = adata[adata.obs[label_key].isin(keep)].copy()
#        adata.obs[label_key] = adata.obs[label_key].cat.reorder_categories(dendrogram, ordered=True)

    # Read NSForest results and align to cluster order
    df = pd.read_csv(results_csv)
    for col in (, "NSForest_markers"):
        if col not in df.columns:
            raise ValueError(f"Column '{col}' missing in {results_csv}")
    df = df[df["cluster_header"].isin(dendrogram)].copy()
    df["cluster_header"l] = df["cluster_header"].astype("category").cat.set_categories(dendrogram)
    df = df.sort_values("cluster_header")

    # Optionally trim markers per cluster
    if top_n is not None and top_n > 0:
        def _trim(v: str) -> str:
            if not isinstance(v, str):
                return v
            s = v.strip()
            if s.startswith("[") and s.endswith("]"):
                try:
                    lst = [str(x).strip() for x in ast.literal_eval(s)]
                except Exception:
                    lst = [p.strip() for p in s.strip("[]").split(",")]
            else:
                lst = [p.strip() for p in s.split(";")]
            return ";".join(lst[:top_n])
        df["NSForest_markers"] = df["NSForest_markers"].map(_trim)

    # Build markers_dict from results
    def _split(val: str) -> list[str]:
        if not isinstance(val, str):
            return []
        s = val.strip()
        if s.startswith("[") and s.endswith("]"):
            try:
                lst = list(ast.literal_eval(s))
                return [str(x).strip() for x in lst if str(x).strip()]
            except Exception:
                pass
        return [g.strip() for g in s.split(";") if g.strip()]

    markers_dict: Dict[str, List[str]] = {
        row["cluster_header"]: _split(row["NSForest_markers"]) for _, row in df.iterrows() if row["cluster_header"] in dendrogram
    }

    # Union of plotted genes in dendrogram cluster order (stable, no duplicates)
    union_order: List[str] = []
    seen = set()
    for cl in dendrogram:
        for g in markers_dict.get(cl, []):
            if g not in seen:
                seen.add(g)
                union_order.append(g)
    if not union_order:
        raise ValueError("No marker genes found to plot (empty union from results).")

    # Keep only genes present in the AnnData
    present = [g for g in union_order if g in adata.var_names]
    if not present:
        raise ValueError("None of the requested marker genes are present in adata.var_names.")

    # Slice to present genes; ensure plotting uses X, not .raw; apply log1p if requested
    ad_for_plot = adata[:, present].copy()
    ad_for_plot.raw = None
    if log1p:
        X = ad_for_plot.X
        if sparse.issparse(X):
            ad_for_plot.X = sparse.csr_matrix(np.log1p(X.toarray()))
        else:
            ad_for_plot.X = np.log1p(X)

    # Render via NSForest plotting (no internal saving) and grab the figure
    ax = ns.pl.stackedviolin(
        ad_for_plot,
        markers_dict,
        label_key,
        dendrogram=dendrogram,
        save=False,
        output_folder=".",
        outputfilename_suffix=".",
    )
    fig = ax.get_figure() if hasattr(ax, "get_figure") else plt.gcf()

    # Save exactly what was requested; no dpi override
    if png_out:
        Path(png_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(png_out), bbox_inches="tight")
    if svg_out:
        Path(svg_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(svg_out), bbox_inches="tight", format="svg")

    plt.close(fig)
    return None

