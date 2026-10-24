"""
Compare NS-Forest results with a published results CSV (for example the Jorstad et al. 2023 files on Zenodo).

For each cluster of either file it records whether the markers are the same (as a set and in order)
and the difference in each metric. Exact agreement is expected when the input data, the NS-Forest
version and the library versions are the same, so every difference is listed.

Saves:
  comparison_{label}.csv      one row per cluster
  comparison_{label}.txt      the summary
"""
import ast
import pandas as pd

from .common_utils import log_section, logger

COUNT_COLUMNS = ["clusterSize", "TN", "FP", "FN", "TP"]
METRIC_TOLERANCE = 1e-6


def _markers_column(df):
    for name in ("NSForest_markers", "markers"):
        if name in df.columns:
            return name
    raise ValueError(f"No markers column in {list(df.columns)}")


def _load(path):
    df = pd.read_csv(path)
    col = _markers_column(df)
    df["_markers"] = [ast.literal_eval(m) if isinstance(m, str) else list(m) for m in df[col]]
    if "PPV" not in df.columns and "precision" in df.columns:
        df["PPV"] = df["precision"]
    # recall from the counts, so a file that lacks the column (the published subclass CSV) compares too
    df["_recall"] = df["TP"] / (df["TP"] + df["FN"])
    return df.set_index("clusterName")


def run_compare(ours_csv, published_csv, label):
    """Compare two results files cluster by cluster; the summary is also logged."""
    log_section(f"Compare with published results: {label}")
    ours, pub = _load(ours_csv), _load(published_csv)

    rows = []
    for name in sorted(set(ours.index) | set(pub.index)):
        row = {"clusterName": name, "in_ours": name in ours.index, "in_published": name in pub.index}
        if row["in_ours"] and row["in_published"]:
            o, p = ours.loc[name], pub.loc[name]
            row["markers_ours"] = ", ".join(o["_markers"])
            row["markers_published"] = ", ".join(p["_markers"])
            row["markers_same_set"] = set(o["_markers"]) == set(p["_markers"])
            row["markers_same_order"] = list(o["_markers"]) == list(p["_markers"])
            for col in ("f_score", "PPV"):
                row[f"{col}_diff"] = float(o[col]) - float(p[col])
            row["recall_diff"] = float(o["_recall"]) - float(p["_recall"])
            for col in COUNT_COLUMNS:
                if col in ours.columns and col in pub.columns:
                    row[f"{col}_diff"] = int(o[col]) - int(p[col])
        rows.append(row)
    report = pd.DataFrame(rows)

    both = report[report["in_ours"] & report["in_published"]]
    diff_cols = [c for c in report.columns if c.endswith("_diff")]
    metric_ok = (both[[c for c in diff_cols if c.split("_diff")[0] in ("f_score", "PPV", "recall")]].abs() <= METRIC_TOLERANCE).all(axis=1)
    count_ok = (both[[c for c in diff_cols if c.split("_diff")[0] in COUNT_COLUMNS]] == 0).all(axis=1)
    exact = both["markers_same_order"] & metric_ok & count_ok

    lines = [
        f"{label}",
        f"clusters in ours: {int(report['in_ours'].sum())}   in published: {int(report['in_published'].sum())}   in both: {len(both)}",
        f"only in ours: {list(report.loc[~report['in_published'], 'clusterName'])}",
        f"only in published: {list(report.loc[~report['in_ours'], 'clusterName'])}",
        f"same markers (as a set): {int(both['markers_same_set'].sum())} of {len(both)}",
        f"same markers in the same order: {int(both['markers_same_order'].sum())} of {len(both)}",
        f"exact (markers, counts, metrics within {METRIC_TOLERANCE}): {int(exact.sum())} of {len(both)}",
        f"largest f_score difference: {both['f_score_diff'].abs().max() if len(both) else 'n/a'}",
        f"clusters that are not exact: {list(both.loc[~exact, 'clusterName'])}",
    ]
    summary = "\n".join(lines)
    for line in lines:
        logger.info(line)

    report.to_csv(f"comparison_{label}.csv", index=False)
    with open(f"comparison_{label}.txt", "w") as f:
        f.write(summary + "\n")
    logger.info(f"Saved: comparison_{label}.csv, comparison_{label}.txt")
