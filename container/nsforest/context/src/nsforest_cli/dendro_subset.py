# nsforest_cli/dendro_subset.py
from __future__ import annotations
from pathlib import Path
from typing import List, Optional

import scanpy as sc
import pandas as pd
import numpy as np

def _ordered_categories_from_uns(uns: dict) -> List[str]:
    # Try common keys that Scanpy stores
    for k in ("categories_ordered", "ordered_categories"):
        if k in uns and uns[k] is not None:
            return list(uns[k])
    # Fallback: derive from categories + index order
    if "categories" in uns and "categories_idx_ordered" in uns:
        cats = list(map(str, uns["categories"]))
        idx  = list(map(int, uns["categories_idx_ordered"]))
        return [cats[i] for i in idx]
    # Last resort: try dendrogram_info["ivl"] (labels produced by scipy)
    info = uns.get("dendrogram_info", {})
    if "ivl" in info and info["ivl"] is not None:
        return list(map(str, info["ivl"]))
    raise ValueError("Could not determine ordered categories from .uns dendrogram payload.")

def subset_by_dendro_range_run(
    h5ad_in: Path,
    h5ad_out: Path,
    label_key: str,
    start: int,
    end: int,
    invert: bool = False,
) -> None:
    """
    Subset cells by taking clusters whose positions in the dendrogram leaf order
    fall within [start, end] (0-based, inclusive). If invert=True, exclude that range.
    """
    adata = sc.read_h5ad(str(h5ad_in))

    key = f"dendrogram_{label_key}"
    if key not in adata.uns:
        # Compute dendrogram if missing (uses PCA fallback automatically)
        sc.tl.dendrogram(adata, groupby=label_key)

    order = _ordered_categories_from_uns(adata.uns[key])

    n = len(order)
    if start < 0 or end < 0 or start >= n or end >= n:
        raise ValueError(f"Range out of bounds: start={start}, end={end}, n_leaves={n}")
    if end < start:
        raise ValueError(f"end ({end}) must be >= start ({start})")

    selected = set(order[start:end + 1])

    col = pd.Series(adata.obs[label_key].astype("string"))
    keep = col.isin(selected)
    if invert:
        keep = ~keep

    adata_sub = adata[keep.to_numpy()].copy()
    # Ensure clean categories on the grouping column
    adata_sub.obs[label_key] = adata_sub.obs[label_key].astype("category")
    adata_sub.obs[label_key] = adata_sub.obs[label_key].cat.remove_unused_categories()

    adata_sub.write_h5ad(str(h5ad_out))


def print_dendro_order_run(h5ad_in: Path, label_key: str) -> List[str]:
    """
    Return the dendrogram leaf order (list of label strings).
    """
    adata = sc.read_h5ad(str(h5ad_in))
    key = f"dendrogram_{label_key}"
    if key not in adata.uns:
        sc.tl.dendrogram(adata, groupby=label_key)
    return _ordered_categories_from_uns(adata.uns[key])

