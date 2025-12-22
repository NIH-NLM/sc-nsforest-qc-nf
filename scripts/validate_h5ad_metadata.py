#!/usr/bin/env python

import pandas as pd
import scanpy as sc
import argparse
import os

def validate_h5ad_metadata(csv_path, log_path):
    df = pd.read_csv(csv_path)

    def log(msg):
        print(msg)
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(msg + "\n")

    for idx, row in df.iterrows():
        log(f"\n=== {row['author']} ({row['publication']}, {row['publication_date']}) ===")
        h5ad_path = row['h5ad_file']
        log(f"Loading: {h5ad_path}")

        if h5ad_path.startswith("http"):
            log("Skipped: Remote .h5ad file detected (not downloaded)")
            continue

        if not os.path.exists(h5ad_path):
            log("File not found")
            continue

        try:
            adata = sc.read_h5ad(h5ad_path)
        except Exception as e:
            log(f"Failed to load: {e}")
            continue

        log(f"Shape: {adata.shape}")
        log(str(adata))

        log("\n.obs keys:")
        log(str(adata.obs.columns.tolist()))

        log("\n.obsm keys:")
        log(str(list(adata.obsm.keys())))

        label_key = row['label_key']
        embedding_key = row['embedding_key']

        if label_key in adata.obs:
            log(f"Label key '{label_key}' found with {adata.obs[label_key].nunique()} unique values")
        else:
            log(f"Label key '{label_key}' NOT found")

        if embedding_key in adata.obsm:
            log(f"Embedding key '{embedding_key}' found")
        else:
            log(f"Embedding key '{embedding_key}' NOT found")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="CSV with h5ad_file, label_key, and embedding_key columns")
    parser.add_argument("output", help="Log file to write output to")
    args = parser.parse_args()

    validate_h5ad_metadata(args.csv, args.output)

