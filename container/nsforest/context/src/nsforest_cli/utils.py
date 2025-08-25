import ast
from pathlib import Path
import pandas as pd
from typing import List, Dict

def parse_markers_field(x: str) -> List[str]:
    """Parse NSForest_markers field into a list of ENSG IDs."""
    if x is None:
        return []
    s = str(x).strip()
    if not s:
        return []
    if s.startswith("[") and s.endswith("]"):
        try:
            lst = ast.literal_eval(s)
            return [str(v).strip().strip("'\"") for v in lst if str(v).strip()]
        except Exception:
            pass
    s = s.strip("[]")
    parts = []
    for token in s.replace("'", "").replace('"', "").split(","):
        token = token.strip()
        if token:
            parts.extend(token.split())
    return [p for p in parts if p]

def markers_from_nsforest_results(
    results_csv: Path,
    cluster_col: str = "clusterName",
    markers_col: str = "NSForest_markers",
) -> Dict[str, List[str]]:
    """Read results.csv and return {cluster -> [ENSG,...]} preserving order."""
    df = pd.read_csv(results_csv, dtype=str)
    if cluster_col not in df.columns or markers_col not in df.columns:
        raise ValueError(f"Expected columns '{cluster_col}' and '{markers_col}' in {results_csv}")

    markers: Dict[str, List[str]] = {}
    for _, row in df.iterrows():
        cl = str(row[cluster_col])
        raw = row[markers_col]
        lst = parse_markers_field(raw)
        if lst:
            markers[cl] = lst
    return markers

