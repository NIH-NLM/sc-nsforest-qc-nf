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
    log_section("NSForest: Run NSForest")
    
    output_folder = create_output_dir(organ, first_author, year) + "/"
    outputfilename_prefix = cluster_header.replace(" ", "_")
    
    # Load prepared adata (positive genes filtered by prep_medians_binary_scores)
    adata_prep = load_h5ad(h5ad_path, cluster_header)

    # Restore varm keys from h5py-safe names back to original names.
    # prep_medians_binary_scores.py replaced '/' with '__SLASH__' before writing.
    # nsforesting.NSForest() requires the original keys e.g. 'medians_subclass.full'.
    _SLASH = "__SLASH__"
    medians_key       = 'medians_' + cluster_header
    binary_scores_key = 'binary_scores_' + cluster_header
    safe_medians_key       = medians_key.replace("/", _SLASH)
    safe_binary_scores_key = binary_scores_key.replace("/", _SLASH)

    if safe_medians_key in adata_prep.varm:
        adata_prep.varm[medians_key]       = adata_prep.varm[safe_medians_key]
        del adata_prep.varm[safe_medians_key]
    if safe_binary_scores_key in adata_prep.varm:
        adata_prep.varm[binary_scores_key] = adata_prep.varm[safe_binary_scores_key]
        del adata_prep.varm[safe_binary_scores_key]

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
        save=False,
        save_supplementary=False,
        output_folder=output_folder,
        outputfilename_prefix=outputfilename_prefix
    )
    
    logger.info(f"NSForest results shape: {results.shape}")
    
    output_csv = f"{output_folder}{outputfilename_prefix}_results.csv"
    results.to_csv(output_csv, index=False)
    logger.info(f"Saved: {output_csv}")
    logger.info("NSForest complete!")
