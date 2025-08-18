from pathlib import Path
import json
import pandas as pd
import scanpy as sc
import nsforest as ns

def eval_markers_run(h5ad_in: Path, markers_json: Path, label_key: str) -> pd.DataFrame:
    """
    Evaluate user marker sets using nsforest.evaluating helpers.
    Returns a DataFrame; caller writes CSV.
    """
    adata = sc.read_h5ad(str(h5ad_in))
    markers = json.loads(markers_json.read_text())  # {cluster: [genes], ...} or similar
    prepared = ns.utils.prepare_markers(markers) if hasattr(ns, "utils") else markers
    df = ns.evaluating.DecisionTree(adata, label_key, prepared)
    return df

