# container/nsforest/context/src/nsforest_cli/dotplot.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import scanpy as sc
import matplotlib.pyplot as plt


def _slice_leaves(
    leaves: list[str],
    *,
    leaf_range: Optional[str] = None,
    leaf_indices: Optional[Sequence[int]] = None,
) -> list[str]:
    """Return a subset of leaves by indices or by slice string 'start:end'.
    - If both are provided, `leaf_indices` wins.
    - `leaf_range` uses Python-style slicing semantics (end exclusive).
    """
    if leaf_indices:
        idx = [i for i in leaf_indices if -len(leaves) <= i < len(leaves)]
        return [leaves[i] for i in idx]
    if leaf_range:
        try:
            start_s, end_s = leaf_range.split(":", 1)
            start = int(start_s) if start_s.strip() != "" else None
            end = int(end_s) if end_s.strip() != "" else None
            return leaves[slice(start, end)]
        except Exception:
            # fall back to full order on parse issues
            return leaves
    return leaves


def dotplot_run(
    h5ad_in: Path,
    results_csv: Path,
    *,
    label_key: str,
    cluster_col: str = "clusterName",
    markers_col: str = "NSForest_markers",
    clusters: Optional[Sequence[str]] = None,
    top_n: Optional[int] = None,
    png_out: Optional[Path] = None,
    svg_out: Optional[Path] = None,
    log1p: bool = True,
    use_ensembl: bool = True,
    id_source: str = "var_names",
    # accept subsetting hints for parity with main.py
    leaf_range: Optional[str] = None,
    leaf_indices: Optional[Sequence[int]] = None,
):
    """Matrix/dot plot using nsforest's plotting API, with dendrogram order.

    - No directory creation; save exactly the filenames requested.
    - Honors `log1p` by transforming a copy of AnnData before plotting.
    - Supports optional dendrogram subsetting via `leaf_indices` / `leaf_range`
      and/or explicit `clusters` list (intersection applied).
    """
    import pandas as pd
    import ast
    import numpy as np
    from scipy import sparse
    import nsforest as ns

    adata = sc.read_h5ad(str(h5ad_in))
    if label_key not in adata.obs:
        raise ValueError(f"obs['{label_key}'] not found in {h5ad_in}")
    adata.obs[label_key] = adata.obs[label_key].astype("category")

    # Ensure dendrogram exists; get leaves
    dendro_key = f"dendrogram_{label_key}"
    if dendro_key not in adata.uns or "categories_ordered" not in adata.uns.get(dendro_key, {}):
        sc.tl.dendrogram(adata, groupby=label_key)
    leaves = list(adata.uns[dendro_key]["categories_ordered"])

    # Apply leaf slicing and/or explicit clusters
    sliced = _slice_leaves(leaves, leaf_range=leaf_range, leaf_indices=leaf_indices)
    if clusters:
        want = set(clusters)
        dendrogram = [c for c in sliced if c in want]
    else:
        dendrogram = sliced

    # Subset AnnData to selected categories (if not full)
    if set(dendrogram) != set(leaves):
        adata = adata[adata.obs[label_key].isin(dendrogram)].copy()
    # Enforce ordered categories for plotting
    adata.obs[label_key] = adata.obs[label_key].cat.set_categories(dendrogram, ordered=True)

    # Load results and align to dendrogram
    df = pd.read_csv(results_csv)
    for col in (cluster_col, markers_col):
        if col not in df.columns:
            raise ValueError(f"Column '{col}' missing in {results_csv}")
    df = df[df[cluster_col].isin(dendrogram)].copy()
    df[cluster_col] = df[cluster_col].astype("category").cat.set_categories(dendrogram)
    df = df.sort_values(cluster_col)

    # Optional top-N trim
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
        df[markers_col] = df[markers_col].map(_trim)

    # Build markers_dict
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
        row[cluster_col]: _split(row[markers_col]) for _, row in df.iterrows() if row[cluster_col] in dendrogram
    }

    # Prepare plotting copy and apply log1p if requested
    ad_for_plot = adata.copy()
    if log1p:
        X = ad_for_plot.X
        if sparse.issparse(X):
            ad_for_plot.X = sparse.csr_matrix(np.log1p(X.toarray()))
        else:
            ad_for_plot.X = np.log1p(X)

    # Render via nsforest plotting (NO saving here) to capture a figure
    ax = ns.pl.matrixplot(
        ad_for_plot,
        markers_dict,
        label_key,
        dendrogram=dendrogram,
        save=False,
        output_folder=".",
        outputfilename_suffix=".",
    )
    fig = ax.get_figure() if hasattr(ax, "get_figure") else plt.gcf()

    # Save
    if png_out:
        Path(png_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(png_out), bbox_inches="tight")
    if svg_out:
        Path(svg_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(svg_out), bbox_inches="tight", format="svg")

    plt.close(fig)
    return None

