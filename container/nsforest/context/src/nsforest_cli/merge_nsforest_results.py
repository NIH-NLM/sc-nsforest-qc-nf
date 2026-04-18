"""
Merge partial NSForest results files and save csv + pkl.
Also generates supplementary marker files from the merged results.

Corresponds to DEMO_NS-Forest_workflow.py: Section 3 (gather phase)

Saves:
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_results.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_results_symbols.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_results.pkl
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_results_symbols.pkl
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_markers.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_markers_symbols.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_markers_onTarget.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_markers_onTarget_symbols.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_markers_onTarget_supp.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_markers_onTarget_supp_symbols.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_gene_selection.csv
  {organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}_gene_selection_symbols.csv
"""

import ast
import os

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

    # --- Merge symbol-keyed partials ---
    sym_dfs = []
    for filepath in partial_files:
        sym_path = str(filepath).replace('.csv', '_symbols.csv')
        if not os.path.exists(sym_path):
            continue
        try:
            df = pd.read_csv(sym_path)
            if df.empty:
                logger.warning(f"Skipping empty symbol partial: {sym_path}")
                continue
            sym_dfs.append(df)
        except pd.errors.EmptyDataError:
            logger.warning(f"Skipping empty symbol partial: {sym_path}")
            continue
        
    if sym_dfs:
        results_symbols = pd.concat(sym_dfs, axis=0, ignore_index=True)
        results_symbols.to_csv(f"{prefix}_results_symbols.csv", index=False)
        results_symbols.to_pickle(f"{prefix}_results_symbols.pkl")
        logger.info(f"Saved: {prefix}_results_symbols.csv")
        logger.info(f"Saved: {prefix}_results_symbols.pkl")
    else:
        logger.warning("No _symbols partial files found — skipping symbol merge")
        results_symbols = None

    # --- Generate supplementary marker files ---
    def _write_marker_files(results_df, suffix=""):
        """
        Emit markers.csv, markers_onTarget.csv, markers_onTarget_supp.csv, gene_selection.csv
        for a given results DataFrame. `suffix` is inserted before `.csv`.
        """

        markers_df = results_df[['clusterName', 'NSForest_markers', 'f_score']].copy()
        markers_df.to_csv(f"{prefix}_markers{suffix}.csv", index=False)
        logger.info(f"Saved: {prefix}_markers{suffix}.csv")

        if 'onTarget' in results_df.columns:
            ontarget_df = results_df[results_df['onTarget'] > 0][
                ['clusterName', 'NSForest_markers', 'f_score', 'onTarget', 'precision', 'recall']
            ].copy()
            ontarget_df.to_csv(f"{prefix}_markers_onTarget{suffix}.csv", index=False)
            logger.info(f"Saved: {prefix}_markers_onTarget{suffix}.csv")

            ontarget_supp = results_df[results_df['onTarget'] > 0].copy()
            ontarget_supp.to_csv(f"{prefix}_markers_onTarget_supp{suffix}.csv", index=False)
            logger.info(f"Saved: {prefix}_markers_onTarget_supp{suffix}.csv")

        all_markers = []
        for _, row in results_df.iterrows():
            markers = row['NSForest_markers']
            if pd.isna(markers):
                continue
            if isinstance(markers, str):
                markers = ast.literal_eval(markers)

            for gene in markers:
                all_markers.append({'clusterName': row['clusterName'], 'gene': gene})

                gene_sel_df = pd.DataFrame(all_markers)
                gene_sel_df.to_csv(f"{prefix}_gene_selection{suffix}.csv", index=False)
                logger.info(f"Saved: {prefix}_gene_selection{suffix}.csv")

            # ENSG outputs (unchanged behavior)
            _write_marker_files(results, "")

            # Symbol outputs (new)
            if results_symbols is not None:
                _write_marker_files(results_symbols, "_symbols")

            logger.info("Merge NSForest results complete!")

