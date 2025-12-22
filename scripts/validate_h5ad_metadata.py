#!/usr/bin/env python

import argparse
import pandas as pd
import scanpy as sc
import nsforest as ns
import os
import urllib.request
import tempfile

def validate_h5ad_metadata(csv_path, output_log_path):
    df = pd.read_csv(csv_path)

    with open(output_log_path, 'w') as f_log:
        for _, row in df.iterrows():
            try:
                author = row['author']
                h5ad_url = row['h5ad']
                cluster_header = row['cluster_header']
                label_key = row.get('label_key', None)

                f_log.write(f"\n=== {author} ===\n")
                f_log.write(f"Loading: {h5ad_url}\n")

                if h5ad_url.startswith("http"):
                    with tempfile.NamedTemporaryFile(suffix=".h5ad") as tmp_file:
                        urllib.request.urlretrieve(h5ad_url, tmp_file.name)
                        adata = sc.read_h5ad(tmp_file.name)
                else:
                    adata = sc.read_h5ad(h5ad_url)

                f_log.write(str(adata) + "\n")

                if label_key and label_key in adata.obs.columns:
                    f_log.write(f"Overriding cluster_header with label_key: {label_key}\n")
                    cluster_header = label_key

                if cluster_header not in adata.obs.columns:
                    f_log.write(f"Cluster header '{cluster_header}' not found in obs.\n")
                    continue

                n_clusters = len(adata.obs[cluster_header].unique())
                fig_width = int(n_clusters / 5)
                fig_height = max([2, int(max([len(z) for z in adata.obs[cluster_header].unique()]) / 30) + 1])

                output_folder = os.path.dirname(h5ad_url) if not h5ad_url.startswith("http") else "."
                if not os.path.exists(output_folder):
                    os.makedirs(output_folder)

                h5ad_basename = os.path.basename(h5ad_url)
                outputfilename_prefix = h5ad_basename.replace(".h5ad", "")
                outputfilename_suffix = ""

                ns.pp.dendrogram(
                    adata,
                    cluster_header,
                    figsize=(fig_width, fig_height),
                    tl_kwargs={"optimal_ordering": True},
                    save="svg",
                    output_folder=output_folder,
                    outputfilename_suffix=outputfilename_suffix
                )

                df_cluster_sizes = pd.DataFrame(adata.obs[cluster_header].value_counts())
                df_cluster_sizes.to_csv(os.path.join(output_folder, f"{outputfilename_prefix}_cluster_sizes.csv"))

                cluster_order = [x.strip() for x in adata.uns["dendrogram_" + cluster_header]['categories_ordered']]
                pd.DataFrame({'cluster_order': cluster_order}).to_csv(
                    os.path.join(output_folder, f"{outputfilename_prefix}_cluster_order.csv"), index=False
                )

                df_normal = pd.DataFrame({
                    'n_obs': [adata.n_obs],
                    'n_vars': [adata.n_vars],
                    'n_clusters': [n_clusters]
                })
                df_normal.to_csv(os.path.join(output_folder, f"{outputfilename_prefix}_summary_normal.csv"), index=False)

                f_log.write("Metadata validation and dendrogram completed.\n")

            except Exception as e:
                f_log.write(f"Error processing {author}: {str(e)}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="Input CSV with metadata")
    parser.add_argument("output", help="Output log file path")
    args = parser.parse_args()
    validate_h5ad_metadata(args.csv, args.output)

if __name__ == "__main__":
    main()

