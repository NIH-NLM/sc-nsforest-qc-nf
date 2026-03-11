"""
Compute cluster statistics.

Corresponds to DEMO_NS-forest_workflow.ipynb: Section 2 - Cluster statistics
"""

import numpy as np
import pandas as pd

from .common_utils import (
    get_output_prefix,
    load_h5ad,
    log_section,
    logger,
    save_dataframe,
    setup_file_logging
)


def compute_cluster_statistics(adata, cluster_header):
    """
    Compute basic statistics for each cluster.
    
    Args:
        adata: AnnData object
        cluster_header: Column name for clusters
        
    Returns:
        DataFrame with cluster statistics
    """
    
    logger.info("Computing cluster statistics...")
    
    total_cells = adata.n_obs
    cluster_counts = adata.obs[cluster_header].value_counts().sort_index()
    
    stats_df = pd.DataFrame({
        'cluster': cluster_counts.index,
        'n_cells': cluster_counts.values,
        'percentage': (cluster_counts.values / total_cells * 100).round(2)
    })
    
    logger.info(f"Total cells: {total_cells}")
    logger.info(f"Number of clusters: {len(stats_df)}")
    logger.info(f"Cells per cluster - Min: {stats_df['n_cells'].min()}, "
                f"Max: {stats_df['n_cells'].max()}, "
                f"Mean: {stats_df['n_cells'].mean():.1f}")
    
    return stats_df


def run_cluster_stats(h5ad_path, cluster_header, organ, first_author, year):
    """
    Main function to compute cluster statistics.
    
    Args:
        h5ad_path: Path to h5ad file
        cluster_header: Column name for clusters
        organ: Organ/tissue type
        first_author: First author surname
        year: Publication year
    """

    setup_file_logging("cluster_stats")
    
    log_section("NSForest: Cluster Statistics")
    
    # Create output directory
    output_prefix = get_output_prefix(organ, first_author, year, cluster_header)

    # Load data
    adata = load_h5ad(h5ad_path, cluster_header)
    
    # Compute statistics (exactly as in notebook)
    stats_df = compute_cluster_statistics(adata, cluster_header)
    
    # Save statistics table
    stats_path = f"{output_prefix}_cluster_statistics"
    save_dataframe(stats_df, stats_path, formats=['csv'])
    
    logger.info("Cluster statistics complete!")
