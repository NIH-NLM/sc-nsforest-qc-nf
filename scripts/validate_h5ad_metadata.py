
import os
import sys
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

# Insert local NSForest path
sys.path.insert(0, os.path.abspath("/home/jovyan/session_data/NSForest"))
import nsforest.pp as ns

def validate_h5ad_metadata(csv_path, output_path=None):
    df = pd.read_csv(csv_path)
    outputs = []

    for idx, row in df.iterrows():
        author = row.get('author', f'Entry {idx}')
        h5ad = row.get('h5ad', None)
        label_key = row.get('label_key', None)
        cluster_header = label_key

        if not h5ad or not label_key:
            outputs.append(f"Skipping {author}: missing h5ad or label_key\n")
            continue

        outputs.append(f"\n=== {author} ===")
        outputs.append(f"Loading: {h5ad}")

        try:
            if h5ad.startswith("http"):
                outputs.append(f"Skipped: Remote .h5ad file detected (not downloaded)")
                continue

            adata = sc.read_h5ad(h5ad)
            outputs.append("Successfully loaded h5ad")
            outputs.append(str(adata))

            if cluster_header not in adata.obs:
                outputs.append(f"Missing cluster header: {cluster_header}")
                continue

            output_folder = os.path.dirname(h5ad)
            filename = os.path.basename(h5ad)
            outputfilename_prefix = filename.replace(".h5ad", "")
            outputfilename_suffix = outputfilename_prefix

            # Auto-adjust figsize
            n_clusters = adata.obs[cluster_header].nunique()
            fig_width = int(n_clusters / 5)
            fig_height = max([2, int(max([len(z) for z in adata.obs[cluster_header].unique()]) / 30) + 1])

            outputs.append("Running nsforest.pp.dendrogram...")
            ns.pp.dendrogram(
                adata,
                cluster_header,
                figsize=(fig_width, fig_height),
                tl_kwargs={"optimal_ordering": True},
                save="svg",
                output_folder=output_folder,
                outputfilename_suffix=outputfilename_suffix,
            )

            # Cluster sizes
            df_cluster_sizes = pd.DataFrame(adata.obs[cluster_header].value_counts())
            outputs.append("Cluster sizes:")
            outputs.append(str(df_cluster_sizes))
            df_cluster_sizes.to_csv(os.path.join(output_folder, outputfilename_prefix + "_cluster_sizes.csv"))

            # Cluster order
            cluster_order = [
                x.strip() for x in adata.uns["dendrogram_" + cluster_header]["categories_ordered"]
            ]
            pd.DataFrame({"cluster_order": cluster_order}).to_csv(
                os.path.join(output_folder, outputfilename_prefix + "_cluster_order.csv"), index=False
            )

            # Summary of adata
            df_normal = pd.DataFrame(
                {"n_obs": [adata.n_obs], "n_vars": [adata.n_vars], "n_clusters": [n_clusters]}
            )
            outputs.append("Data summary:")
            outputs.append(str(df_normal))
            df_normal.to_csv(os.path.join(output_folder, outputfilename_prefix + "_summary_normal.csv"), index=False)

        except Exception as e:
            outputs.append(f"Error processing {author}: {str(e)}")

    if output_path:
        with open(output_path, "w") as f:
            f.write("\n".join(outputs))
    else:
        print("\n".join(outputs))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Validate h5ad metadata and extract dendrogram info")
    parser.add_argument("csv", help="CSV file with h5ad info")
    parser.add_argument("output", nargs="?", default=None, help="Optional: write output to file")
    args = parser.parse_args()

    validate_h5ad_metadata(args.csv, args.output)
