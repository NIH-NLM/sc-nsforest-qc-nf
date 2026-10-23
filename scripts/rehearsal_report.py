#!/usr/bin/env python3
"""
Judge a rehearsal run: did every task finish, and are the files that matter there?

  python scripts/rehearsal_report.py <outdir of the rehearsal>

Prints the report and exits 0 when everything that must be there is there and no task failed, otherwise 1.
"""
import csv
import glob
import os
import sys
from collections import Counter

# (label, glob pattern in the output folder, must it be there)
EXPECTED = [
    ("dendrogram figure",                 "dendrogram_*.svg",                              True),
    ("cluster order",                     "cluster_order_*.csv",                           True),
    ("cluster sizes",                     "cluster_sizes_*.csv",                           False),
    ("cluster statistics",                "cluster_statistics*",                           False),
    ("cluster to ontology mapping",       "cluster_cid*",                                  False),
    ("medians (csv)",                     "medians_ensg_*.csv",                            True),
    ("binary scores (csv)",               "binary_scores_ensg_*.csv",                      True),
    ("histograms",                        "hist_nonzero_*.svg",                            False),
    ("subclass NS-Forest results",        "results_ensg_*CrossArea_subclass*.csv",         True),
    ("cluster NS-Forest results",         "results_ensg_*CrossArea_cluster*.csv",          True),
    ("markers",                           "markers_ensg_*.csv",                            False),
    ("gene selection",                    "gene_selection_ensg_*.csv",                     False),
    ("marker combinations",               "results_combinations_*CrossArea_cluster*.csv",  True),
    ("boxplots (svg)",                    "boxplot_*.svg",                                 True),
    ("dotplot (svg)",                     "dotplot_*.svg",                                 True),
    ("stacked violin (svg)",              "stacked_violin_*.svg",                          False),
    ("matrix plot (svg)",                 "matrixplot_*.svg",                              False),
    ("silhouette scores, subclass",       "silhouette_scores_*CrossArea_subclass*.csv",    True),
    ("silhouette scores, cluster",        "silhouette_scores_*CrossArea_cluster*.csv",     True),
    ("silhouette cluster summary",        "scsilhouette_cluster_summary_*.csv",            True),
    ("silhouette annotation",             "annotation_*.json",                             True),
    ("2D projection",                     "*projection*",                                  False),
    ("dataset summary",                   "*dataset_summary*.csv",                         False),
    ("S3 manifest",                       "master_s3_manifest.csv",                        False),
]


def main(outdir):
    ok = True
    lines = []

    trace = os.path.join(outdir, "pipeline_trace.txt")
    if not os.path.exists(trace):
        lines.append(f"NO TRACE FILE: {trace} (the workflow did not get as far as running tasks)")
        ok = False
    else:
        with open(trace) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        statuses = Counter(r["status"] for r in rows)
        lines.append(f"tasks: {len(rows)}   " + "   ".join(f"{k}: {v}" for k, v in sorted(statuses.items())))
        bad = [r for r in rows if r["status"] not in ("COMPLETED", "CACHED")]
        for r in bad:
            lines.append(f"  NOT COMPLETED: {r['name']}  ({r['status']}, exit {r.get('exit', '?')})  work dir hash {r.get('hash', '?')}")
        if bad or not rows:
            ok = False

    lines.append("")
    for label, pattern, must in EXPECTED:
        found = glob.glob(os.path.join(outdir, "**", pattern), recursive=True)
        mark = "ok     " if found else ("MISSING" if must else "absent ")
        lines.append(f"  {mark} {label:32s} {pattern}   ({len(found)} file{'s' if len(found) != 1 else ''})")
        if must and not found:
            ok = False

    lines.append("")
    lines.append("RESULT: PASS" if ok else "RESULT: FAIL  (look at the NOT COMPLETED and MISSING lines; 'absent' lines are files whose names I could not be sure of)")
    print("\n".join(lines))
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))
