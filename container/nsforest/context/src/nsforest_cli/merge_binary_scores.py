"""
Merge partial binary scores files and save csv + pkl.
"""

import pandas as pd
from pathlib import Path

from .common_utils import (
    create_output_dir,
    log_section,
    logger
)


def run_merge_binary_scores(partial_files, cluster_header, organ, first_author, year):
    """
    Merge partial binary scores files and save csv + pkl.
    
    Saves:
    - binary_scores.csv
    - binary_scores.pkl
    """
    log_section("NSForest: Merge Binary Scores")
    
    logger.info(f"Merging {len(partial_files)} partial binary scores files...")
    
    # Read and concatenate
    dfs = []
    for filepath in partial_files:
        df = pd.read_csv(filepath, index_col=0)
        dfs.append(df)
    
    df_binary_scores = pd.concat(dfs, axis=0)
    logger.info(f"Complete binary scores matrix: {df_binary_scores.shape}")
    
    # Save outputs
    output_folder = create_output_dir(organ, first_author, year)
    outputfilename_prefix = cluster_header
    
    df_binary_scores.to_csv(output_folder + "/" + outputfilename_prefix + "_binary_scores.csv")
    df_binary_scores.to_pickle(output_folder + "/" + outputfilename_prefix + "_binary_scores.pkl")
    
    logger.info(f"Saved: {outputfilename_prefix}_binary_scores.csv")
    logger.info(f"Saved: {outputfilename_prefix}_binary_scores.pkl")
    logger.info("Merge binary scores complete!")
