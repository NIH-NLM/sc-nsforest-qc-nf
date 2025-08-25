# container/nsforest/context/src/nsforest_cli/dendro_subset.py
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple
import scanpy as sc
import numpy as np
from anndata import AnnData


def _parse_leaf_range(s: str) -> Tuple[int, int]:
    # accepts "start:end" (python-style, end exclusive)
    s = s.strip()
    if ":" not in s:
        raise ValueError("Leaf range must be 'start:end' with ':'")
    a, b = s.split(":", 1)
    start = int(a) if a else 0
    end = int(b) if b else 10**9
    if end <= start:
        raise ValueError("Leaf range end must be greater than start")
    return start, end


def leaves_from_dendrogram(
    adata: AnnData,
    label_key: str,
    *,
    leaf_range: Optional[str] = None,
    leaf_indices: Optional[Sequence[int]] = None,
) -> List[str]:
    """
    Return cluster labels in dendrogram leaf order, optionally subset by positions.
    """
    # Compute dendrogram (idempotent; overwrites/creates uns)
    sc.tl.dendrogram(adata, groupby=label_key)
    dendro_key = f"dendrogram_{label_key}"
    info = adata.uns.get(dendro_key, {}).get("dendrogram_info", {})
    ivl: List[str] = list(info.get("ivl", []))  # leaf labels in order

    if not ivl:
        raise RuntimeError(f"No dendrogram info found under adata.uns['{dendro_key}']")

    # Full list if no subsetting
    if not leaf_range and not leaf_indices:
        return ivl

    if leaf_range:
        start, end = _parse_leaf_range(leaf_range)
        end = min(end, len(ivl))
        if start < 0 or start >= len(ivl):
            raise IndexError("leaf_range start out of bounds")
        if end <= start:
            raise IndexError("leaf_range end must be > start")
        return ivl[start:end]

    # explicit indices
    idxs = np.array(list(leaf_indices), dtype=int)
    if (idxs < 0).any() or (idxs >= len(ivl)).any():
        raise IndexError("leaf_indices out of bounds for dendrogram leaves")
    return [ivl[i] for i in idxs]

