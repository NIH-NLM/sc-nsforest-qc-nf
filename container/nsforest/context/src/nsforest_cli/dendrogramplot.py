from pathlib import Path
import scanpy as sc
import nsforest as ns

def dendrogramplot_run(h5ad_in: Path, label_key: str):
    """Return AnnData with dendrogram metadata attached via NSForest."""
    adata = sc.read_h5ad(str(h5ad_in)).copy()
    adata = ns.pp.dendrogram(adata, label_key)
    return adata

