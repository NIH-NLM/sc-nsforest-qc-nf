from pathlib import Path
import scanpy as sc
import nsforest as ns

def prep_binary_scores_run(h5ad_in: Path, label_key: str):
    """Return AnnData with binary scores via nsforest.pp.prep_binary_scores."""
    adata = sc.read_h5ad(str(h5ad_in)).copy()
    adata = ns.pp.prep_binary_scores(adata, label_key)
    return adata

