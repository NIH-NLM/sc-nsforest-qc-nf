from pathlib import Path
import scanpy as sc
import nsforest as ns

def prep_medians_run(h5ad_in: Path, label_key: str):
    """Return AnnData with per-cluster medians computed via nsforest.pp.prep_medians."""
    adata = sc.read_h5ad(str(h5ad_in)).copy()
    adata = ns.pp.prep_medians(adata, label_key)
    return adata

