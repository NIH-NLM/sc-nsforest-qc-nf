from pathlib import Path
import re
import scanpy as sc
import pandas as pd

_KEEP = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-")

def _clean_token(s: str) -> str:
    s2 = ''.join(ch if ch in _KEEP else '_' for ch in str(s))
    s2 = re.sub(r'_+', '_', s2).strip('_')
    return s2 or "cluster"

def sanitize_labels_run(h5ad_in: Path, label_key: str):
    """
    Clean labels in .obs[label_key] by replacing non [A-Za-z0-9_-] with '_',
    collapsing repeats, trimming ends. Preserve duplicates across cells.
    If two *distinct* original labels collide to the same cleaned token,
    add a numeric suffix per *original label* (not per cell).
    """
    adata = sc.read_h5ad(str(h5ad_in)).copy()
    if label_key not in adata.obs.columns:
        raise KeyError(f"{label_key!r} not found in .obs")

    # Work on unique original labels to decide final cleaned names
    orig_vals = pd.Index(adata.obs[label_key].astype("string")).unique()

    base_map = {}          # cleaned_base -> count assigned so far
    label_map = {}         # original_label -> final_cleaned

    for orig in orig_vals:
        base = _clean_token(orig)
        if base not in base_map:
            base_map[base] = 0
            final = base
        else:
            base_map[base] += 1
            final = f"{base}_{base_map[base]}"
        label_map[orig] = final

    # Apply the mapping back to all cells; make categorical; drop unused cats
    new_col = pd.Series(adata.obs[label_key].astype("string")).map(label_map)

    # Preserve appearance order of categories, then drop unused
    cats = pd.Index(pd.unique(new_col))
    adata.obs[label_key] = pd.Categorical(new_col, categories=cats, ordered=False)
    adata.obs[label_key] = adata.obs[label_key].cat.remove_unused_categories()

    return adata

