"""
Generate dendrogram and cluster statistics.

Corresponds to DEMO_NS-Forest_workflow.py: Section 2

Saves:
  {organ}_{first_author}_{year}_{cluster_header}_dendrogram.svg
  {organ}_{first_author}_{year}_{cluster_header}_cluster_sizes.csv
  {organ}_{first_author}_{year}_{cluster_header}_cluster_order.csv
  {organ}_{first_author}_{year}_{cluster_header}_summary_normal.csv
"""

import matplotlib
matplotlib.use("Agg")

import pandas as pd
import nsforest as ns

from .common_utils import (
    load_h5ad,
    log_section,
    logger
)


def run_dendrogram(h5ad_path, cluster_header, organ, first_author, year):
    """
    Generate dendrogram and cluster statistics.
    """
    log_section("NSForest: Dendrogram")

    cluster_header_safe = cluster_header.replace(" ", "_")
    prefix = f"{organ}_{first_author}_{year}_{cluster_header_safe}"

    adata = load_h5ad(h5ad_path, cluster_header)

    n_clusters = adata.obs[cluster_header].nunique()
    logger.info(f"Number of clusters: {n_clusters}")

    # Dendrogram — save svg
    logger.info("Creating dendrogram...")
    ns.pp.dendrogram(
        adata,
        cluster_header,
        tl_kwargs={'optimal_ordering': True},
        save="svg",
        output_folder="",
        outputfilename_suffix=prefix
    )

    # Cluster sizes
    df_cluster_sizes = pd.DataFrame(adata.obs[cluster_header].value_counts())
    df_cluster_sizes.to_csv(f"{prefix}_cluster_sizes.csv")
    logger.info(f"Saved: {prefix}_cluster_sizes.csv")

    # Cluster order
    cluster_order = [x.strip() for x in adata.uns["dendrogram_" + cluster_header]['categories_ordered']]
    pd.DataFrame({'cluster_order': cluster_order}).to_csv(f"{prefix}_cluster_order.csv", index=False)
    logger.info(f"Saved: {prefix}_cluster_order.csv")

    # Summary statistics
    df_normal = pd.DataFrame({'n_obs': [adata.n_obs], 'n_vars': [adata.n_vars], 'n_clusters': [n_clusters]})
    df_normal.to_csv(f"{prefix}_summary_normal.csv", index=False)
    logger.info(f"Saved: {prefix}_summary_normal.csv")

    logger.info("Dendrogram complete!")
