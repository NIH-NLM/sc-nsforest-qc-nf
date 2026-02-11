"""
Run NSForest algorithm to identify marker genes (parallelized by cluster).

Uses nsforesting.NSForest() to find minimum marker gene combinations.
Each cluster is processed independently using cluster_list parameter.
"""

import pandas as pd
from nsforest import nsforesting

from .common_utils import (
    create_output_dir,
    load_h5ad,
    log_section,
    logger
)


def run_nsforest(h5ad_path, cluster_header, organ, first_author, year, 
                 cluster_list=None, n_trees=1000, n_genes_eval=6):
    """
    Run NSForest algorithm for specified cluster(s).
    
    Requires adata_prep with:
    - varm['medians_{cluster_header}']
    - varm['binary_scores_{cluster_header}']
    
    When cluster_list has one cluster: saves partial results with unique name.
    When cluster_list is None: processes all clusters, saves complete results.
    
    Parameters
    ----------
    cluster_list : list of str, optional
        Specific clusters to process (for parallelization)
    n_trees : int
        Number of trees in random forest (default: 1000)
    n_genes_eval : int
        Number of top genes to evaluate for combinations (default: 6)
    """
    log_section("NSForest: Run NSForest")
    
    output_folder = create_output_dir(organ, first_author, year) + "/"
    outputfilename_prefix = cluster_header
    
    # Load prepared adata (must have medians and binary_scores in varm)
    adata_prep = load_h5ad(h5ad_path, cluster_header)
    
    # Verify required data
    medians_key = 'medians_' + cluster_header
    binary_scores_key = 'binary_scores_' + cluster_header
    
    if medians_key not in adata_prep.varm:
        raise ValueError(f"Missing {medians_key} in adata.varm. Run prep_medians first.")
    if binary_scores_key not in adata_prep.varm:
        raise ValueError(f"Missing {binary_scores_key} in adata.varm. Run prep_binary_scores first.")
    
    logger.info(f"Medians shape: {adata_prep.varm[medians_key].shape}")
    logger.info(f"Binary scores shape: {adata_prep.varm[binary_scores_key].shape}")
    
    # Run NSForest
    if cluster_list is not None:
        logger.info(f"Running NSForest for cluster(s): {cluster_list}")
    else:
        logger.info("Running NSForest for all clusters")
    
    results = nsforesting.NSForest(
        adata_prep, 
        cluster_header,
        cluster_list=cluster_list if cluster_list else [],
        n_trees=n_trees,
        n_genes_eval=n_genes_eval,
        save=False,  # We'll save manually for parallelization
        save_supplementary=False,
        output_folder=output_folder,
        outputfilename_prefix=outputfilename_prefix
    )
    
    logger.info(f"NSForest results shape: {results.shape}")
    
    # Save results
    if cluster_list is not None and len(cluster_list) == 1:
        # Partial results - unique filename
        cluster_safe = cluster_list[0].replace(' ', '_').replace('/', '-')
        output_csv = f"{output_folder}{outputfilename_prefix}_results_{cluster_safe}.csv"
    else:
        # Complete results
        output_csv = f"{output_folder}{outputfilename_prefix}_results.csv"
    
    results.to_csv(output_csv, index=False)
    logger.info(f"Saved: {output_csv}")
    logger.info("NSForest complete!")
