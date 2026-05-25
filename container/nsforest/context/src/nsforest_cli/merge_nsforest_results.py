"""
Merge partial NSForest results files and save csv + pkl.
Also generates supplementary marker files from the merged results.

Corresponds to DEMO_NS-Forest_workflow.py: Section 3 (gather phase)

Saves:
  results_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.csv
  results_symbols_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.csv
  results_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.pkl
  results_symbols_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.pkl
  markers_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.csv
  markers_symbols_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.csv
  markers_onTarget_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.csv
  markers_onTarget_symbols_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.csv
  markers_onTarget_supp_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.csv
  markers_onTarget_symbols_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.csv
  gene_selection_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.csv
  gene_selection_symbols_{organ}_{first_author}_{year}_{cluster_header}_{embedding}_{vid}.csv
"""

import ast
import csv
import pandas as pd

from .common_utils import (
    get_output_prefix,
    load_h5ad,
    log_section,
    logger
)

def _write_results(results_df, prefix, suffix="", adata=None):
    """Write results csv/pkl and the four marker files for one flavor.

    If ``adata`` is supplied and has an ``adata.var['gene_symbol']`` column, the
    ``NSForest_markers`` column is first converted from ENSG lists to gene-symbol
    lists. Unmapped genes stay as ENSG. If ``adata`` is supplied but the column is
    missing, the symbol-flavored outputs are skipped with a warning — we don't
    silently write ENSG markers into files tagged ``_symbols``.
    """
    if adata is not None:
        if 'gene_symbol' not in adata.var.columns:
            logger.warning(
                f"adata.var['gene_symbol'] missing — skipping '{suffix}' outputs"
            )
            return
        sym_map = dict(zip(adata.var_names, adata.var['gene_symbol']))
        results_df = results_df.copy()
        results_df['NSForest_markers'] = [
            (
                [sym_map.get(g, g) for g in (ast.literal_eval(m) if isinstance(m, str) else m)]
                if not pd.isna(m) else m
            )
            for m in results_df['NSForest_markers']
        ]

    results_df.to_csv(f"results{suffix}_{prefix}.csv", index=False, quoting=csv.QUOTE_NONE)
    results_df.to_pickle(f"results{suffix}_{prefix}.pkl")
    logger.info(f"Saved: results{suffix}_{prefix}.csv")
    logger.info(f"Saved: results{suffix}_{prefix}.pkl")

    markers_df = results_df[['clusterName', 'NSForest_markers', 'f_score']].copy()
    markers_df.to_csv(f"markers{suffix}_{prefix}.csv", index=False, quoting=csv.QUOTE_NONE)
    logger.info(f"Saved: markers{suffix}_{prefix}.csv")

    if 'onTarget' in results_df.columns:
        ontarget_df = results_df[results_df['onTarget'] > 0][
            ['clusterName', 'NSForest_markers', 'f_score', 'onTarget', 'precision', 'recall']
        ].copy()
        ontarget_df.to_csv(f"markers_onTarget{suffix}_{prefix}.csv", index=False, quoting=csv.QUOTE_NONE)
        logger.info(f"Saved: markers_onTarget{suffix}_{prefix}.csv")

        ontarget_supp = results_df[results_df['onTarget'] > 0].copy()
        ontarget_supp.to_csv(f"markers_onTarget_supp{suffix}_{prefix}.csv", index=False, quoting=csv.QUOTE_NONE)
        logger.info(f"Saved: markers_onTarget_supp{suffix}_{prefix}.csv")

    all_markers = []
    for _, row in results_df.iterrows():
        markers = row['NSForest_markers']
        if not isinstance(markers, (list, str)) and pd.isna(markers):
            continue
        if isinstance(markers, str):
            markers = ast.literal_eval(markers)
        for gene in markers:
            all_markers.append({'clusterName': row['clusterName'], 'gene': gene})

    gene_sel_df = pd.DataFrame(all_markers)
    gene_sel_df.to_csv(f"gene_selection{suffix}_{prefix}.csv", index=False, quoting=csv.QUOTE_NONE)
    logger.info(f"Saved: gene_selection{suffix}_{prefix}.csv")

def run_merge_nsforest_results(
        partial_files,
        filtered_h5ad,
        cluster_header,
        organ,
        first_author,
        journal,
        year,
        embedding,
        dataset_version_id,
):
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

     # ENSG flavor
    _write_results(results, prefix)

    # Symbol flavor
    if not results.empty:
        adata = load_h5ad(filtered_h5ad, cluster_header)
        _write_results(results, prefix, suffix="_symbols", adata=adata)

    logger.info(f"Complete results: {results.shape}")
    logger.info("Merge and write NSForest results complete!")

