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

    logger.info("Creating dendrogram...")

    import matplotlib.pyplot as plt
    cluster_names = adata.obs[cluster_header].astype(str).unique()
    n_clusters    = len(cluster_names)
    max_label_len = max((len(c) for c in cluster_names), default=10)

    logger.info(f"Number of clusters: {n_clusters}")

    # Width grows with cluster count; height grows if long labels need rotation room.
    # Tune these coefficients against your worst-case dataset.
    fig_width  = max(10.0, n_clusters * 0.45)
    fig_height = max(6.0, 4.0 + max_label_len * 0.12)
    figsize    = (fig_width, fig_height)
    logger.info(f"Dendrogram figsize: {figsize}")

    with plt.rc_context({
            'figure.figsize':      (fig_width, fig_height),
            'savefig.bbox':        'tight',
            'savefig.pad_inches':  0.4,
    }):
        try:
            ns.pp.dendrogram(
                adata,
                cluster_header,
                tl_kwargs={"optimal_ordering": True},
                pl_kwargs={"show": False},
                save=False,
                figsize=figsize
            )
            plt.savefig(f"{prefix}_dendrogram.svg")
            plt.close()
            
        except ValueError as e:
            if "negative distances" in str(e):
                logger.warning(f"Optimal ordering failed — retrying without: {e}")
                ns.pp.dendrogram(
                    adata,
                    cluster_header,
                    tl_kwargs={"optimal_ordering": False},
                    pl_kwargs={"show": False},
                    save=False,
                    figsize=figsize
                )
                plt.savefig(f"{prefix}_dendrogram.svg")
                plt.close()

            else:
                raise

    # Cluster sizes
    cluster_counts   = adata.obs[cluster_header].value_counts()
    df_cluster_sizes = pd.DataFrame({
        "cluster": cluster_counts.index.tolist(),
        "count":   cluster_counts.values
    })
    df_cluster_sizes.to_csv(f"{prefix}_cluster_sizes.csv", index=False)
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
