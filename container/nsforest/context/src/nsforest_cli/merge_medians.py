"""
Merge partial median files and save csv + pkl.

Corresponds to DEMO_NS-forest_workflow.py: Section 3 (gather phase)
DEMO: df_medians.to_csv() and df_medians.to_pickle()
"""

import pandas as pd
from pathlib import Path

from .common_utils import (
    create_output_dir,
    log_section,
    logger
)


def run_merge_medians(partial_files, cluster_header, organ, first_author, year):
    """
    Merge partial median files and save csv + pkl.
    
    DEMO Section 3:
    df_medians.to_csv(output_folder + outputfilename_prefix + "_medians.csv")
    df_medians.to_pickle(output_folder + outputfilename_prefix + "_medians.pkl")
    """
    log_section("NSForest: Merge Medians")
    
    logger.info(f"Merging {len(partial_files)} partial median files...")
    
    # Read and concatenate
    dfs = []
    for filepath in partial_files:
        df = pd.read_csv(filepath, index_col=0)
        dfs.append(df)
    
    df_medians = pd.concat(dfs, axis=0)
    logger.info(f"Complete median matrix: {df_medians.shape}")
    
    # Save outputs (matching DEMO)
    output_folder = create_output_dir(organ, first_author, year)
    outputfilename_prefix = cluster_header
    
    df_medians.to_csv(output_folder + "/" + outputfilename_prefix + "_medians.csv")
    df_medians.to_pickle(output_folder + "/" + outputfilename_prefix + "_medians.pkl")
    
    logger.info(f"Saved: {outputfilename_prefix}_medians.csv")
    logger.info(f"Saved: {outputfilename_prefix}_medians.pkl")
    logger.info("Merge medians complete!")
