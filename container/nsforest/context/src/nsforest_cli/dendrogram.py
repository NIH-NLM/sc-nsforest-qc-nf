"""
Generate dendrogram and cluster statistics.

Corresponds to DEMO_NS-Forest_workflow.py: Section 2

Saves:
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_dendrogram.svg
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_cluster_sizes.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_cluster_order.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_summary_normal.csv
"""

import matplotlib
matplotlib.use("Agg")

import pandas as pd
import nsforest as ns
import os

from .common_utils import (
    get_output_prefix,
    load_h5ad,
    log_section,
    logger
)


def run_dendrogram(h5ad_path, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id):
    """
    Generate dendrogram and cluster statistics.
    """
    log_section("NSForest: Dendrogram")

    prefix = get_output_prefix( organ, first_author, journal, year, cluster_header, embedding, dataset_version_id )
    logger.info(f"Dendrogram prefix: {prefix}")

    adata = load_h5ad(h5ad_path, cluster_header)

    n_clusters = adata.obs[cluster_header].nunique()
    logger.info(f"Number of clusters: {n_clusters}")

    logger.info("Creating dendrogram...")
    # Scale figure height with number of clusters so rotated labels are not clipped
    fig_width = max(12, n_clusters * 0.5)
    fig_height = max(4, n_clusters * 0.2)
    figsize = (fig_width, fig_height)
    logger.info(f"Dendrogram figsize: {figsize}")
    try:
        ns.pp.dendrogram(
            adata, cluster_header,
            tl_kwargs={"optimal_ordering": True}, save="svg",
            output_folder="", outputfilename_suffix=prefix,
            figsize=figsize
        )
    except ValueError as e:
        if "negative distances" in str(e):
            logger.warning(f"Optimal ordering failed — retrying without: {e}")
            ns.pp.dendrogram(
                adata, cluster_header,
                tl_kwargs={"optimal_ordering": False}, save="svg",
                output_folder="", outputfilename_suffix=prefix,
                figsize=figsize
            )
        else:
            raise

    # NSForest library saves as _{prefix}.svg — rename to {prefix}_dendrogram.svg
    src = f"_{prefix}.svg"
    dst = f"{prefix}_dendrogram.svg"
    if os.path.exists(src):
        os.rename(src, dst)
        logger.info(f"Renamed: {src} -> {dst}")
    else:
        logger.warning(f"Expected dendrogram SVG not found: {src}")


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
