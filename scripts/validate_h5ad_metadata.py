#!/usr/bin/env python

import pandas as pd
import scanpy as sc
import argparse
import os

def log(msg, output_log):
    print(msg)
    with open(output_log, "a", encoding="utf-8") as f:
        f.write(msg + "\n")

def validate_h5ad_metadata(row, output_log):
    author = row['author']
    pub = row['publication']
    date = row['publication_date']
    h5ad_path = row['h5ad_file']
    label_key = row['label_key']
    embedding_key = row['embedding_key']

    log(f"\n=== {author} ({pub}, {date}) ===", output_log)
    log(f"Loading: {h5ad_path}", output_log)

    if h5ad_path.startswith("http"):
        log("Skipped: Remote .h5ad file detected (not downloaded)", output_log)
        return

    if not os.path.exists(h5ad_path):
        log("File not found", output_log)
        return

    try:
        adata = sc.read_h5ad(h5ad_path)
    except Exception as e:
        log(f"Failed to load: {e}", output_log)
        return

    log(f"Shape: {adata.shape}", output_log)
    log(str(adata), output_log)

    log("\n.obs keys:", output_log)
    log(str(adata.obs.columns.tolist()), output_log)

    log("\n.obsm keys:", output_log)
    log(str(list(adata.obsm.keys())), output_log)

    if label_key not in adata.obs:
        log(f"Label key '{label_key}' NOT found in adata.obs", output_log)
        return

    if embedding_key not in adata.obsm:
        log(f"Embedding key '{embedding_key}' NOT found in adata.obsm", output_log)

    cluster_header = label_key
    n_clusters = adata.obs[cluster_header].nunique()
    log(f"Assigned cluster_header = '{cluster_header}'", output_log)
    log(f"Number of unique clusters: {n_clusters}", output_log)

    output_folder = os.path.join(os.path.dirname(h5ad_path), "nsforest_output")
    os.makedirs(output_folder, exist_ok=True)
    h5ad_basename = os.path.basename(h5ad_path).replace(".h5ad", "")
    outputfilename_prefix = os.path.join(output_folder, f"{h5ad_basename}_{cluster_header}")
    outputfilename_suffix = "dendrogram"

    fig_width = int(n_clusters / 5)
    fig_height = max([2, int(max([len(z) for z in adata.obs[cluster_header].unique()]) / 30) + 1])

    # Dendrogram
    pp.dendrogram(
        adata,
        cluster_header,
        figsize=(fig_width, fig_height),
        tl_kwargs={"optimal_ordering": True},
        save="svg",
        output_folder=output_folder,
        outputfilename_suffix=outputfilename_suffix
    )

    # Cluster sizes
    df_cluster_sizes = pd.DataFrame(adata.obs[cluster_header].value_counts())
    df_cluster_sizes.to_csv(outputfilename_prefix + "_cluster_sizes.csv")
    log("Saved: cluster_sizes.csv", output_log)

    # Cluster order
    cluster_order = [x.strip() for x in adata.uns["dendrogram_" + cluster_header]['categories_ordered']]
    pd.DataFrame({"cluster_order": cluster_order}).to_csv(outputfilename_prefix + "_cluster_order.csv", index=False)
    log("Saved: cluster_order.csv", output_log)

    # Summary
    df_normal = pd.DataFrame({
        "n_obs": [adata.n_obs],
        "n_vars": [adata.n_vars],
        "n_clusters": [n_clusters]
    })
    df_normal.to_csv(outputfilename_prefix + "_summary_normal.csv", index=False)
    log("Saved: summary_normal.csv", output_log)

def main():
    parser = argparse.ArgumentParser(
        description="Inspect AnnData metadata and pre-NS-Forest inputs"
    )
    parser.add_argument("csv",    help="CSV file with h5ad_file, label_key, embedding_key")
    parser.add_argument("output", help="Output log file")
    args = parser.parse_args()

    validate_h5ad_metadata(args.csv, args.output)

    df = pd.read_csv(args.csv_file)
    if os.path.exists(args.output):
        os.remove(args.output)

    for _, row in df.iterrows():
        validate_h5ad_metadata(row, args.output)

if __name__ == "__main__":
    main()


