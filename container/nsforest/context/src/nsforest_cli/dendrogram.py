"""
Generate dendrogram and cluster statistics using NSForest.

"""

import pandas as pd
import nsforest as ns

from .common_utils import (
    create_output_dir,
    load_h5ad,
    log_section,
    logger
)


def run_dendrogram(h5ad_path, cluster_header, organ, first_author, year):
    """
    
    Creates:
    1. dendrogram.svg - hierarchical clustering visualization
    2. cluster_sizes.csv - cell counts per cluster
    3. cluster_order.csv - ordered cluster list from dendrogram
    4. summary_normal.csv - dataset statistics (n_obs, n_vars, n_clusters)
    
    All outputs required for downstream FRMatch analysis.
    """
    log_section("NSForest: Section 2 - Clusters")
    
    # Match DEMO variable names exactly
    output_folder = create_output_dir(organ, first_author, year) + "/"
    outputfilename_suffix = cluster_header
    outputfilename_prefix = cluster_header
    
    adata = load_h5ad(h5ad_path, cluster_header)
    
    # Number of clusters 
    n_clusters = adata.obs[cluster_header].nunique()
    logger.info(f"Number of clusters: {n_clusters}")
    
    # Auto-adjust figsize 
    fig_width = int(n_clusters / 5)
    fig_height = max([2, int(max([len(z) for z in adata.obs[cluster_header].unique()]) / 30) + 1])
    
    # Dendrogram and save svg 
    logger.info("Creating dendrogram...")
    ns.pp.dendrogram(
        adata, 
        cluster_header, 
        figsize=(fig_width, fig_height),
        tl_kwargs={'optimal_ordering': True},
        save="svg",
        output_folder=output_folder,
        outputfilename_suffix=outputfilename_suffix
    )
    logger.info(f"Saved: {outputfilename_suffix}_dendrogram.svg")
    
    # Cluster sizes 
    df_cluster_sizes = pd.DataFrame(adata.obs[cluster_header].value_counts())
    df_cluster_sizes.to_csv(output_folder + outputfilename_prefix + "_cluster_sizes.csv")
    logger.info(f"Saved: {outputfilename_prefix}_cluster_sizes.csv")
    
    # Cluster order 
    cluster_order = [x.strip() for x in adata.uns["dendrogram_" + cluster_header]['categories_ordered']]
    pd.DataFrame({'cluster_order': cluster_order}).to_csv(
        output_folder + outputfilename_prefix + "_cluster_order.csv", 
        index=False
    )
    logger.info(f"Saved: {outputfilename_prefix}_cluster_order.csv")
    
    # Summary statistics 
    df_normal = pd.DataFrame({
        'n_obs': [adata.n_obs], 
        'n_vars': [adata.n_vars], 
        'n_clusters': [n_clusters]
    })
    df_normal.to_csv(output_folder + outputfilename_prefix + "_summary_normal.csv", index=False)
    logger.info(f"Saved: {outputfilename_prefix}_summary_normal.csv")
    
    logger.info("Section 2 complete!")
