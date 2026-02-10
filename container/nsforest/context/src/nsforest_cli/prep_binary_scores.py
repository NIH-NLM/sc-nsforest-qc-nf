"""
Compute binary scores for marker discovery.

Binary scores quantify how well each gene distinguishes a cluster from all others.
Uses Mann-Whitney U test or t-test to compare expression in target cluster vs rest.

Corresponds to DEMO_NS-forest_workflow.ipynb: Section 4 - Binary scores
"""

import numpy as np
import pandas as pd
from scipy import stats
from typing import Literal

from .common_utils import (
    create_output_dir,
    get_output_prefix,
    log_section,
    logger
)


def compute_binary_scores(median_df, method='mwu'):
    """
    Compute binary scores for each cluster.
    
    For each cluster, compares median expression in that cluster vs all other clusters
    for every gene. Higher scores indicate better cluster-specific markers.
    
    Parameters
    ----------
    median_df : pandas.DataFrame
        Complete median expression matrix (clusters × genes)
    method : {'mwu', 'ttest'}
        Statistical test method:
        - 'mwu': Mann-Whitney U test (non-parametric, default)
        - 'ttest': Student's t-test (parametric)
        
    Returns
    -------
    pandas.DataFrame
        Binary scores matrix (clusters × genes)
        Higher values = better markers for that cluster
        
    Notes
    -----
    Binary score calculation:
    1. For each cluster C and gene G:
       - Group 1: Expression in cluster C
       - Group 2: Expression in all other clusters
    2. Compute test statistic (U or t)
    3. Convert to score (higher = more specific to cluster)
    """
    logger.info(f"Computing binary scores using {method.upper()} test...")
    logger.info(f"Median matrix shape: {median_df.shape}")
    
    clusters = median_df.index.tolist()
    genes = median_df.columns.tolist()
    
    n_clusters = len(clusters)
    n_genes = len(genes)
    
    logger.info(f"Processing {n_clusters} clusters × {n_genes} genes...")
    
    # Initialize scores matrix
    binary_scores = np.zeros((n_clusters, n_genes))
    
    for i, cluster in enumerate(clusters):
        if (i + 1) % 10 == 0 or i == 0:
            logger.info(f"  Processing cluster {i+1}/{n_clusters}: {cluster}")
        
        # Get expression for this cluster and all others
        cluster_expr = median_df.loc[cluster].values
        other_expr = median_df.drop(cluster).values
        
        for j, gene in enumerate(genes):
            gene_in_cluster = cluster_expr[j]
            gene_in_others = other_expr[:, j]
            
            if method == 'mwu':
                # Mann-Whitney U test
                # Higher U when cluster has higher expression
                try:
                    u_stat, p_val = stats.mannwhitneyu(
                        [gene_in_cluster], 
                        gene_in_others,
                        alternative='greater'
                    )
                    binary_scores[i, j] = u_stat
                except:
                    binary_scores[i, j] = 0
                    
            elif method == 'ttest':
                # T-test
                try:
                    t_stat, p_val = stats.ttest_ind(
                        [gene_in_cluster],
                        gene_in_others,
                        equal_var=False
                    )
                    binary_scores[i, j] = t_stat if t_stat > 0 else 0
                except:
                    binary_scores[i, j] = 0
    
    # Convert to DataFrame
    scores_df = pd.DataFrame(
        binary_scores,
        index=clusters,
        columns=genes
    )
    
    logger.info(f"Binary scores computed: {scores_df.shape}")
    logger.info(f"Score range: [{scores_df.min().min():.2f}, {scores_df.max().max():.2f}]")
    
    return scores_df


def run_prep_binary_scores(median_matrix_path, cluster_header, organ, first_author, year, method='mwu'):
    """
    Main function to compute binary scores.
    
    Parameters
    ----------
    median_matrix_path : str or Path
        Path to complete median expression matrix CSV
    cluster_header : str
        Column name for clusters
    organ : str
        Organ/tissue type
    first_author : str
        First author surname
    year : str
        Publication year
    method : {'mwu', 'ttest'}
        Statistical test method (default: 'mwu')
        
    Output Files
    ------------
    Creates: outputs_{organ}_{first_author}_{year}/{cluster_header}_binary_scores.csv
    Format: CSV with clusters as rows, genes as columns, binary scores as values
    
    Examples
    --------
    >>> run_prep_binary_scores(
    ...     median_matrix_path='outputs_kidney_Lake_2023/subclass.full_medians.csv',
    ...     cluster_header='subclass.full',
    ...     organ='kidney',
    ...     first_author='Lake',
    ...     year='2023',
    ...     method='mwu'
    ... )
    """
    log_section("NSForest: Binary Scores")
    
    logger.info(f"Method: {method.upper()}")
    logger.info(f"Loading median matrix: {median_matrix_path}")
    
    # Load median matrix
    median_df = pd.read_csv(median_matrix_path, index_col=0)
    logger.info(f"Median matrix shape: {median_df.shape} (clusters × genes)")
    
    # Compute binary scores
    scores_df = compute_binary_scores(median_df, method=method)
    
    # Save binary scores
    output_dir = create_output_dir(organ, first_author, year)
    output_prefix = get_output_prefix(output_dir, cluster_header)
    output_file = f"{output_prefix}_binary_scores.csv"
    
    scores_df.to_csv(output_file)
    logger.info(f"Saved binary scores: {output_file}")
    logger.info("Binary scores complete!")
    logger.info("")
    logger.info("Next step: Plot histograms (Section 5)")
