# nsforest_cli/dendrogramplot_run.py
from __future__ import annotations

from pathlib import Path
from typing import Optional, List

import matplotlib.pyplot as plt
import scanpy as sc
import nsforest as ns
import typer
from nsforest import ns, nsforesting, utils, NSFOREST_VERSION

# We expect `ensembl_lookup(ensg_id) -> str | list | dict | None`.
try:
    from nsforest_cli.ensembl_lookup import ensembl_lookup  # type: ignore
except Exception:  # pragma: no cover
    ensembl_lookup = None  # type: ignore

def dendrogramplot_run(
        *,
        h5ad_in:       Path,
        label_key:     str,
        h5ad_out:      Path,
        png_out:       Optional[Path],
        svg_out:       Optional[Path],
        leaf_range:    Optional[str],
        leaf_indices:  Optional[List[int]],):
        
    """Populate `uns[dendrogram_<label_key>]` and set `var_names` via `ensembl_lookup`.

    Rules (tutorial-aligned):
    - For each ENSG in original `var_names`, call `ensembl_lookup(ENSG)`:
      - if lookup resolves to exactly one gene symbol → use that as the new var name;
      - otherwise → keep the original ENSG ID.
    - Mutates AnnData in-place; optional `h5ad_out` persists the result.
    """
    if ensembl_lookup is None:  # pragma: no cover
        typer.echo("ensembl_lookup function not available; ensure it is importable.", err=True)
        raise typer.Exit(code=2)

    adata = sc.read_h5ad(str(h5ad_in))

    if label_key not in adata.obs:
        typer.echo(f"obs['{label_key}'] not found in AnnData.", err=True)
        raise typer.Exit(code=2)

    # 1) Resolve display names with ensembl_lookup, fallback to ENSG if not uniquely resolved
    def _lookup_unique_symbol(ensg: str) -> Optional[str]:
        try:
            res = ensembl_lookup(ensg)
        except Exception:
            return None
        if res is None:
            return None
        if isinstance(res, str):
            s = res.strip()
            return s or None
        if isinstance(res, (list, tuple, set)):
            # Unique only if it collapses to a single non-empty symbol
            vals = [str(x).strip() for x in res if str(x).strip()]
            uniq = sorted(set(vals))
            return uniq[0] if len(uniq) == 1 else None
        if isinstance(res, dict):
            cand = res.get("gene_symbol") or res.get("symbol") or res.get("name")
            if isinstance(cand, str):
                s = cand.strip()
                return s or None
            if isinstance(cand, (list, tuple, set)):
                vals = [str(x).strip() for x in cand if str(x).strip()]
                uniq = sorted(set(vals))
                return uniq[0] if len(uniq) == 1 else None
        return None

    original = list(adata.var_names)
    new_names = [(_lookup_unique_symbol(ensg) or ensg) for ensg in original]
    adata.var_names = new_names

    # 2) Persist dendrogram structure in `uns` (no files written here)
    ns.pl.dendrogram(adata, cluster_header=label_key, save=True)

    # because we are using this as part of a workflow - we save the persisted dendrogram in a new h5ad
    adata.write_h5ad(str(h5ad_out))

    return None

