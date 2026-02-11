"""
Compute binary scores per cluster (parallelized by cluster).

Binary scores measure how specifically a gene is expressed in one cluster vs others.
Uses ns.pp.prep_binary_scores() on filtered adata from prep_medians.
"""

import pandas as pd
import nsforest as ns

from .common_utils import (
    create_output_dir,
    load_h5ad,
    log_section,
    logger
)

def run_prep_binary_scores(h5ad_path, cluster_header, organ, first_author, year, cluster_list=None):
    """
    Compute binary scores for specified cluster(s).
    
    Saves both partial binary_scores CSV AND adata_prep.h5ad (with binary scores added to varm).
    """
    log_section("NSForest: Prep Binary Scores")
    
    output_folder = create_output_dir(organ, first_author, year)
    outputfilename_prefix = cluster_header
    
    # Load adata_prep (should already be filtered by prep_medians)
    adata_prep = load_h5ad(h5ad_path, cluster_header)
    
    # Run NSForest prep_binary_scores
    logger.info("Running ns.pp.prep_binary_scores()...")
    adata_prep = ns.pp.prep_binary_scores(adata_prep, cluster_header)
    
    # Save adata_prep (all parallel jobs create identical adata_prep)
    adata_prep_path = f"{output_folder}/adata_prep.h5ad"
    adata_prep.write_h5ad(adata_prep_path)
    logger.info(f"Saved: adata_prep.h5ad")
    
    # Extract binary scores from varm
    df_binary_scores = adata_prep.varm['binary_scores_' + cluster_header].T
    
    # Filter to specific cluster(s) if requested
    if cluster_list is not None:
        logger.info(f"Filtering to cluster(s): {cluster_list}")
        df_binary_scores = df_binary_scores.loc[cluster_list]
    
    logger.info(f"Binary scores shape: {df_binary_scores.shape}")
    
    # Save with unique filename if single cluster
    if cluster_list is not None and len(cluster_list) == 1:
        cluster_safe = cluster_list[0].replace(' ', '_').replace('/', '-')
        output_csv = f"{output_folder}/{outputfilename_prefix}_binary_scores_{cluster_safe}.csv"
    else:
        output_csv = f"{output_folder}/{outputfilename_prefix}_binary_scores.csv"
    
    df_binary_scores.to_csv(output_csv)
    logger.info(f"Saved: {output_csv}")
    logger.info("Prep binary scores complete!")

