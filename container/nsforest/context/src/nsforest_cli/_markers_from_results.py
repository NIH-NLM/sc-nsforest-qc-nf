from pathlib import Path
from typing import Dict, List
import ast
import pandas as pd

def markers_from_nsforest_results(
    results_csv: Path,
    cluster_col: str = "clusterName",
    markers_col: str = "NSForest_markers",  # or "binary_genes"
) -> Dict[str, List[str]]:
    """
    Build {cluster: [genes]} from NS-Forest results.csv that already contains list-like
    marker columns (e.g., 'NSForest_markers' or 'binary_genes').
    """
    df = pd.read_csv(results_csv)
    if cluster_col not in df.columns:
        raise ValueError(f"'{cluster_col}' not found in results.csv")
    if markers_col not in df.columns:
        raise ValueError(f"'{markers_col}' not found in results.csv")

    out: Dict[str, List[str]] = {}
    for _, row in df[[cluster_col, markers_col]].iterrows():
        cluster = str(row[cluster_col])
        raw = row[markers_col]
        # rows look like "['ENSG...', 'ENSG...']" -> parse safely
        if isinstance(raw, str) and raw.startswith("["):
            genes = [str(g) for g in ast.literal_eval(raw)]
        elif isinstance(raw, list):
            genes = [str(g) for g in raw]
        else:
            # fallback: treat as a single marker token
            genes = [str(raw)]
        out[cluster] = genes
    return out

