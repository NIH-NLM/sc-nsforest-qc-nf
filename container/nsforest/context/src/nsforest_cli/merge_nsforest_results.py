"""
Merge partial NSForest results files and save csv + pkl.
"""

import pandas as pd
from pathlib import Path

from .common_utils import (
    create_output_dir,
    log_section,
    logger
)


def run_merge_nsforest_results(partial_files, cluster_header, organ, first_author, year):
    """
    Merge partial NSForest results files and save csv + pkl.
    
    Saves:
    - results.csv
    - results.pkl
    """
    log_section("NSForest: Merge NSForest Results")
    
    logger.info(f"Merging {len(partial_files)} partial NSForest results files...")
    
    # Read and concatenate
    dfs = []
    for filepath in partial_files:
        df = pd.read_csv(filepath)
        dfs.append(df)
    
    results = pd.concat(dfs, axis=0, ignore_index=True)
    logger.info(f"Complete results: {results.shape} rows")
    
    # Save outputs
    output_folder = create_output_dir(organ, first_author, year)
    outputfilename_prefix = cluster_header
    
    results.to_csv(output_folder + "/" + outputfilename_prefix + "_results.csv", index=False)
    results.to_pickle(output_folder + "/" + outputfilename_prefix + "_results.pkl")
    
    logger.info(f"Saved: {outputfilename_prefix}_results.csv")
    logger.info(f"Saved: {outputfilename_prefix}_results.pkl")
    logger.info("Merge NSForest results complete!")
