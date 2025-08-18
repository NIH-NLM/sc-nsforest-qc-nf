from pathlib import Path
from typing import List, Literal
import pandas as pd
import scanpy as sc

Mode = Literal["exact", "contains", "regex"]
NApol = Literal["drop", "keep", "match"]

def filter_by_obs_run(
    h5ad_in: Path,
    obs_key: str,
    values: List[str],
    *,
    mode: Mode = "exact",          # exact | contains | regex
    case_insensitive: bool = True,
    na_policy: NApol = "drop",     # drop | keep | match
    invert: bool = False,          # keep non-matching rows if True
):
    """
    Return a new AnnData filtered by a *single* .obs[obs_key].
    - values: one or more strings to match
    - mode:   'exact' equality, 'contains' substring, or 'regex'
    - case_insensitive: lowercases both sides if True
    - na_policy: how to treat NA in obs_key before matching
    - invert: keep rows that do NOT match
    """
    if not values:
        raise ValueError("Provide at least one value to match.")
    adata = sc.read_h5ad(str(h5ad_in))
    if obs_key not in adata.obs.columns:
        raise KeyError(f"obs_key '{obs_key}' not in adata.obs")

    s = adata.obs[obs_key].astype("string")
    if na_policy == "drop":
        mask_valid = s.notna()
        s = s[mask_valid]
    elif na_policy == "keep":
        mask_valid = pd.Series(True, index=s.index)
    elif na_policy == "match":
        s = s.fillna("NA")
        mask_valid = pd.Series(True, index=s.index)
    else:
        raise ValueError("na_policy must be 'drop' | 'keep' | 'match'")

    if case_insensitive:
        s_cmp = s.str.lower()
        vals = [v.lower() for v in values]
    else:
        s_cmp = s
        vals = values

    if mode == "exact":
        mask_match = pd.Series(False, index=s_cmp.index)
        for v in vals:
            mask_match |= (s_cmp == v)
    elif mode == "contains":
        mask_match = pd.Series(False, index=s_cmp.index)
        for v in vals:
            mask_match |= s_cmp.str.contains(v, regex=False)
    elif mode == "regex":
        pattern = "(" + "|".join(vals) + ")"
        mask_match = s_cmp.str.contains(pattern, regex=True)
    else:
        raise ValueError("mode must be 'exact' | 'contains' | 'regex'")

    final_mask = mask_valid & ( ~mask_match if invert else mask_match )
    return adata[final_mask].copy()

