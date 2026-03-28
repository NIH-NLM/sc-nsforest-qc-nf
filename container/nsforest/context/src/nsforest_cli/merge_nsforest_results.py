"""
Merge partial NSForest results files and save csv + pkl.
Also generates supplementary marker files from the merged results.

Corresponds to DEMO_NS-Forest_workflow.py: Section 3 (gather phase)

Saves:
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_results.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_results.pkl
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_markers.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_markers_onTarget.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_markers_onTarget_supp.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_gene_selection.csv
"""

import ast

import pandas as pd

from .common_utils import (
    get_output_prefix,
    log_section,
    logger
)


def run_merge_nsforest_results(partial_files, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id):
    """
    Merge partial NSForest results CSV files and save csv + pkl + marker files.
    """
    log_section("NSForest: Merge NSForest Results")

    prefix = get_output_prefix( organ, first_author, journal, year, cluster_header, embedding, dataset_version_id )

    logger.info(f"Merging {len(partial_files)} partial NSForest results files...")

    dfs = []
    for filepath in partial_files:
        try:
            df = pd.read_csv(filepath)
            if df.empty:
                logger.warning(f"Skipping empty partial file: {filepath}")
                continue
            dfs.append(df)
        except pd.errors.EmptyDataError:
            logger.warning(f"Skipping empty partial file: {filepath}")
            continue

    if not dfs:
        logger.warning("No non-empty partial files found — writing empty results")
        results = pd.DataFrame()
    else:
        results = pd.concat(dfs, axis=0, ignore_index=True)

    logger.info(f"Complete results: {results.shape}")

    results.to_csv(f"{prefix}_results.csv", index=False)
    results.to_pickle(f"{prefix}_results.pkl")
    logger.info(f"Saved: {prefix}_results.csv")
    logger.info(f"Saved: {prefix}_results.pkl")

    # --- Generate supplementary marker files ---

    # markers.csv — all clusters with marker genes and scores
    markers_df = results[['clusterName', 'NSForest_markers', 'f_score']].copy()
    markers_df.to_csv(f"{prefix}_markers.csv", index=False)
    logger.info(f"Saved: {prefix}_markers.csv")

    # markers_onTarget.csv — clusters with onTarget > 0
    if 'onTarget' in results.columns:
        ontarget_df = results[results['onTarget'] > 0][
            ['clusterName', 'NSForest_markers', 'f_score', 'onTarget', 'precision', 'recall']
        ].copy()
        ontarget_df.to_csv(f"{prefix}_markers_onTarget.csv", index=False)
        logger.info(f"Saved: {prefix}_markers_onTarget.csv")

        # markers_onTarget_supp.csv — all columns for on-target clusters
        ontarget_supp = results[results['onTarget'] > 0].copy()
        ontarget_supp.to_csv(f"{prefix}_markers_onTarget_supp.csv", index=False)
        logger.info(f"Saved: {prefix}_markers_onTarget_supp.csv")

    # gene_selection.csv — per-cluster gene selection
    all_markers = []
    for _, row in results.iterrows():
        markers = row['NSForest_markers']
        if pd.isna(markers):
            continue
        if isinstance(markers, str):
            markers = ast.literal_eval(markers)
        for gene in markers:
            all_markers.append({'clusterName': row['clusterName'], 'gene': gene})
    gene_sel_df = pd.DataFrame(all_markers)
    gene_sel_df.to_csv(f"{prefix}_gene_selection.csv", index=False)
    logger.info(f"Saved: {prefix}_gene_selection.csv")

    logger.info("Merge NSForest results complete!")
