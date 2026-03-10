"""
Merge partial NSForest results files and save csv + pkl.

Corresponds to DEMO_NS-Forest_workflow.py: Section 3 (gather phase)

Saves:
  {organ}_{first_author}_{year}_{cluster_header}_results.csv
  {organ}_{first_author}_{year}_{cluster_header}_results.pkl
"""

import pandas as pd

from .common_utils import (
    log_section,
    logger
)


def run_merge_nsforest_results(partial_files, cluster_header, organ, first_author, year):
    """
    Merge partial NSForest results CSV files and save csv + pkl.
    """
    log_section("NSForest: Merge NSForest Results")

    cluster_header_safe = cluster_header.replace(" ", "_")
    prefix = f"{organ}_{first_author}_{year}_{cluster_header_safe}"

    logger.info(f"Merging {len(partial_files)} partial NSForest results files...")

    dfs = []
    for filepath in partial_files:
        df = pd.read_csv(filepath)
        dfs.append(df)

    results = pd.concat(dfs, axis=0, ignore_index=True)
    logger.info(f"Complete results: {results.shape}")

    results.to_csv(f"{prefix}_results.csv", index=False)
    results.to_pickle(f"{prefix}_results.pkl")
    logger.info(f"Saved: {prefix}_results.csv")
    logger.info(f"Saved: {prefix}_results.pkl")

    logger.info("Merge NSForest results complete!")
