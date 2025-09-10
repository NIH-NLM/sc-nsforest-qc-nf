# nsforest_cli/symbolize_genes.py

from __future__ import annotations
from pathlib import Path
from typing import Optional
import pandas as pd
import scanpy as sc

def load_mapping(csv_file: Path) -> dict[str, str]:
    df = pd.read_csv(csv_file)
    return dict(zip(df["ensg"], df["symbol"]))

def symbolize_genes_run(
    *,
    h5ad_in: Path,
    h5ad_out: Path,
    symbol_map_csv: Path,
    csv_out: Optional[Path] = None,
):
    """
    Use a local ENSG→symbol CSV to update adata.var_names.
    Stores original ENSGs in .var['orig_var_names'].
    """
    adata = sc.read_h5ad(str(h5ad_in)).copy()

    full_ensgs = list(map(str, adata.var_names))
    stripped_ensgs = [e.split(".")[0] for e in full_ensgs]

    mapping = load_mapping(symbol_map_csv)

    hits = sum(e in mapping for e in stripped_ensgs)
    print(f"[info] Matched {hits}/{len(full_ensgs)} ENSG IDs to symbols")

    adata.var["orig_var_names"] = adata.var_names
    new_var_names = [mapping.get(e, orig) for e, orig in zip(stripped_ensgs, full_ensgs)]

    unresolved = sum(1 for s, o in zip(new_var_names, full_ensgs) if s == o)
    if unresolved > 0:
        print(f"[warn] {unresolved} symbols could not be resolved and remain as ENSG IDs")

    adata.var_names = new_var_names

    if csv_out:
        pd.DataFrame({
            "ensg": full_ensgs,
            "stripped": stripped_ensgs,
            "symbol": new_var_names,
        }).to_csv(csv_out, index=False)

    adata.write_h5ad(str(h5ad_out))

