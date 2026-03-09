"""
Compute median expression per cluster (parallelized by cluster).

Corresponds to DEMO_NS-forest_workflow.py: Section 3 prep
Uses ns.pp.prep_medians() to filter positive genes and compute medians.
"""

import pandas as pd
import nsforest as ns

from .common_utils import (
    create_output_dir,
    get_output_prefix,
    load_h5ad,
    log_section,
    logger
)

def run_prep_medians(h5ad_path, cluster_header, organ, first_author, year, cluster_list=None):
    """
    Compute medians for specified cluster(s).
    
    Saves both partial medians CSV AND adata_prep.h5ad (with positive gene filter applied).
    """
    log_section("NSForest: Prep Medians")
    
    output_folder = create_output_dir(organ, first_author, year)
    outputfilename_prefix = cluster_header.replace(" ", "_")
    
    # Load and prepare data
    adata = load_h5ad(h5ad_path, cluster_header)
    
    # Make a copy
    adata_prep = adata.copy()
    
    # Run NSForest prep_medians (filters positive genes, computes medians)
    logger.info("Running ns.pp.prep_medians()...")
    adata_prep = ns.pp.prep_medians(adata_prep, cluster_header)
    
    # Extract median matrix from varm before writing h5ad
    # (cluster names as column values are content, not h5py keys)
    df_medians = adata_prep.varm['medians_' + cluster_header].T
    del adata_prep.varm['medians_' + cluster_header]

    # Save adata_prep (positive gene filter applied, medians removed from varm)
    adata_prep_path = f"{output_folder}/adata_prep.h5ad"
    adata_prep.write_h5ad(adata_prep_path)
    logger.info(f"Saved: adata_prep.h5ad")
    
    # Filter to specific cluster(s) if requested
    if cluster_list is not None:
        logger.info(f"Filtering to cluster(s): {cluster_list}")
        df_medians = df_medians.loc[cluster_list]
    
    logger.info(f"Median matrix shape: {df_medians.shape}")
    
    output_csv = f"{output_folder}/{outputfilename_prefix}_medians.csv"
    df_medians.to_csv(output_csv)
    logger.info(f"Saved: {output_csv}")
    logger.info("Prep medians complete!")
