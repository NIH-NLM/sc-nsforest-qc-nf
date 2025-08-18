from pathlib import Path
import json
import scanpy as sc
import nsforest as ns

def dotplot_run(h5ad_in: Path, markers_json: Path, label_key: str):
    """Return Matplotlib figure from NSForest dotplot helper."""
    adata = sc.read_h5ad(str(h5ad_in))
    markers = json.loads(Path(markers_json).read_text())
    ax = ns.plotting.dotplot(adata, markers, label_key)
    return ax.figure

