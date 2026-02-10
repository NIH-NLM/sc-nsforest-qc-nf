"""
Prepare median expression matrix per cluster.

This module computes median gene expression for clusters. It supports
parallelization by processing individual clusters or cluster subsets.

Corresponds to DEMO_NS-forest_workflow.ipynb: Section 3 - Prepare medians
"""

import numpy as np
import pandas as pd

from .common_utils import (
    create_output_dir,
    get_output_prefix,
    load_h5ad,
    log_section,
    logger
)


def compute_medians(adata, cluster_header, cluster_list=None):
    """
    Compute median expression for specified clusters.
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with expression data
    cluster_header : str
        Column name in adata.obs containing cluster labels
    cluster_list : list of str, optional
        Specific clusters to compute medians for. If None, compute for all clusters.
        This enables Nextflow parallelization: pass single cluster for per-cluster jobs.
        
    Returns
    -------
    pandas.DataFrame
        Median expression matrix with clusters as rows, genes as columns
        
    Notes
    -----
    For Nextflow parallelization:
    - cluster_list=None: Compute all clusters (sequential, outputs medians.csv)
    - cluster_list=['ClusterA']: Compute single cluster (parallel, outputs medians_partial.csv)
    - Nextflow handles scatter (split clusters) and gather (merge partial CSVs)
    """
    # Determine which clusters to process
    if cluster_list is None:
        clusters_to_process = adata.obs[cluster_header].unique()
        logger.info(f"Computing medians for all {len(clusters_to_process)} clusters")
    else:
        clusters_to_process = cluster_list
        logger.info(f"Computing medians for {len(clusters_to_process)} cluster(s): {clusters_to_process}")
    
    median_dict = {}
    
    for cluster in clusters_to_process:
        logger.info(f"Processing cluster: {cluster}")
        
        # Get cells for this cluster
        cluster_mask = adata.obs[cluster_header] == cluster
        cluster_cells = adata[cluster_mask]
        n_cells = cluster_cells.n_obs
        
        logger.info(f"  Cluster '{cluster}': {n_cells} cells")
        
        # Compute median expression
        if hasattr(cluster_cells, 'X') and cluster_cells.X is not None:
            # Use X matrix (preprocessed data)
            if hasattr(cluster_cells.X, 'toarray'):
                # Sparse matrix
                median_expr = np.median(cluster_cells.X.toarray(), axis=0)
            else:
                # Dense matrix
                median_expr = np.median(cluster_cells.X, axis=0)
        else:
            # Fallback to raw counts
            logger.warning(f"  No X matrix found for cluster '{cluster}', using raw counts")
            if hasattr(cluster_cells.raw.X, 'toarray'):
                median_expr = np.median(cluster_cells.raw.X.toarray(), axis=0)
            else:
                median_expr = np.median(cluster_cells.raw.X, axis=0)
        
        # Ensure 1D array
        median_expr = np.asarray(median_expr).flatten()
        median_dict[cluster] = median_expr
    
    # Create DataFrame: clusters as rows, genes as columns
    median_df = pd.DataFrame(median_dict, index=adata.var_names).T
    
    logger.info(f"Median matrix shape: {median_df.shape} (clusters × genes)")
    logger.info(f"Clusters in matrix: {list(median_df.index)}")
    
    return median_df


def run_prep_medians(h5ad_path, cluster_header, organ, first_author, year, cluster_list=None):
    """
    Main function to prepare median expression matrix.
    
    Parameters
    ----------
    h5ad_path : str or Path
        Path to h5ad file
    cluster_header : str
        Column name for clusters in adata.obs
    organ : str
        Organ/tissue type (e.g., 'kidney', 'heart')
    first_author : str
        First author surname
    year : str
        Publication year
    cluster_list : list of str, optional
        Specific clusters to process. If None, processes all clusters.
        Used by Nextflow for parallelization.
    """
    log_section("NSForest: Prepare Medians")
    
    # Create output directory
    output_dir = create_output_dir(organ, first_author, year)
    output_prefix = get_output_prefix(output_dir, cluster_header)
    
    # Load data
    adata = load_h5ad(h5ad_path, cluster_header)
    
    # Compute median expression
    median_df = compute_medians(adata, cluster_header, cluster_list)

    # Save median matrix with UNIQUE filename per cluster
    if cluster_list is not None and len(cluster_list) == 1:
        # Single cluster - add cluster name to filename for uniqueness
        cluster_safe = cluster_list[0].replace(' ', '_').replace('/', '-')
        output_file = f"{output_prefix}_medians_{cluster_safe}.csv"
        logger.info(f"Note: Partial file for cluster: {cluster_list[0]}")
    else:
        # Multiple clusters or all clusters - use standard name
        output_file = f"{output_prefix}_medians.csv"
    
    median_df.to_csv(output_file)
    logger.info(f"Saved median matrix: {output_file}")
    logger.info("Median preparation complete!")
