from pathlib import Path
import scanpy as sc
import nsforest as ns

def dendrogramplot_run(
    h5ad_in: Path,
    results_csv: Path,   # accepted for CLI symmetry; not used
    label_key: str,
):
    """
    Attach dendrogram metadata. Works whether ns.pp.dendrogram returns AnnData
    or modifies in place and returns None. Falls back to scanpy.tl.dendrogram.
    """
    adata = sc.read_h5ad(str(h5ad_in)).copy()

    # First try NS-Forest helper
    try:
        ret = ns.pp.dendrogram(adata, label_key)
        if ret is not None:
            adata = ret  # some versions return a new AnnData
        # else: modified in place; keep 'adata'
    except Exception:
        # Fallback to scanpy API
        sc.tl.dendrogram(adata, groupby=label_key)

    return adata

