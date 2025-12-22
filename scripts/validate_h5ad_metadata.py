import pandas as pd
import scanpy as sc
import argparse
import os

def validate_h5ad_metadata(csv_path):
    df = pd.read_csv(csv_path)

    for idx, row in df.iterrows():
        print(f"\n=== {row['author']} ({row['publication']}, {row['publication_date']}) ===")
        h5ad_path = row['h5ad_file']
        print(f"Loading: {h5ad_path}")

        # Check for local file
        if h5ad_path.startswith("http"):
            print("Skipped: Remote .h5ad file detected (not downloaded)")
            continue
        if not os.path.exists(h5ad_path):
            print("File not found")
            continue

        # 1. Data 
        try:
            adata = sc.read_h5ad(h5ad_path)
        except Exception as e:
            print(f"Failed to load: {e}")
            continue

        print(f"Shape: {adata.shape}")

        print( adata )

        # 2. Clusters
        #
        # number of clusters
        # n_clusters = adata.obs[cluster_header].nunique()
        # print ( n_clusters )
        
        print("\n.obs keys:")
        print(adata.obs.columns.tolist())

        print("\n.obsm keys:")
        print(list(adata.obsm.keys()))

        label_key = row['label_key']
        embedding_key = row['embedding_key']

        if label_key in adata.obs:
            print(f"Label key '{label_key}' found with {adata.obs[label_key].nunique()} unique values")
        else:
            print(f"Label key '{label_key}' NOT found")

        if embedding_key in adata.obsm:
            print(f"Embedding key '{embedding_key}' found")
        else:
            print(f"Embedding key '{embedding_key}' NOT found")

import sys
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("--output", help="Optional: write output to file")
    args = parser.parse_args()

    if args.output:
        sys.stdout = open(args.output, "w")


