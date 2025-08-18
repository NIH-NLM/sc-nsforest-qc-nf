from pathlib import Path
import re
import scanpy as sc

_KEEP = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-")

def _clean_token(s: str) -> str:
    s2 = ''.join(ch if ch in _KEEP else '_' for ch in s)
    s2 = re.sub(r'_+', '_', s2).strip('_')
    return s2 or "cluster"

def sanitize_labels_run(h5ad_in: Path, label_key: str):
    """
    Collision-safe label sanitizer: replace non [A-Za-z0-9_-] with '_',
    collapse repeats, strip ends, dedupe with suffixes. Returns AnnData.
    """
    adata = sc.read_h5ad(str(h5ad_in)).copy()
    if label_key not in adata.obs.columns:
        raise KeyError(f"label_key '{label_key}' not in adata.obs")
    old = adata.obs[label_key].astype(str).map(_clean_token)

    seen = {}
    uniq = []
    for val in old:
        if val not in seen:
            seen[val] = 0
            uniq.append(val)
        else:
            seen[val] += 1
            uniq.append(f"{val}_{seen[val]}")
    adata.obs[label_key] = uniq
    return adata

