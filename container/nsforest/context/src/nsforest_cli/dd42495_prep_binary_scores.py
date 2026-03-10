from pathlib import Path
import scanpy as sc
import nsforest as ns

def prep_binary_scores_run(
        *,
        h5ad_in: Path,
        label_key: str,
        h5ad_out: Path,
):
    """Return AnnData with binary scores via nsforest.pp.prep_binary_scores."""
    adata = sc.read_h5ad(str(h5ad_in)).copy()
    adata = ns.pp.prep_binary_scores(adata, label_key)

    # because we are using this as part of a workflow - we save the persisted dendrogram in a new h5ad
    adata.write_h5ad(str(h5ad_out))
    
    return None

