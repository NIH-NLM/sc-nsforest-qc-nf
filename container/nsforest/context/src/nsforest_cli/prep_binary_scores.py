"""
Compute binary scores per cluster.

Corresponds to DEMO_NS-Forest_workflow.py: Section 3 prep
Uses ns.pp.prep_medians() then ns.pp.prep_binary_scores() in memory.

Saves:
  {organ}_{first_author}_{year}_{cluster_header}_binary_scores.csv
  {organ}_{first_author}_{year}_{cluster_header}_binary_scores.pkl
"""

import nsforest as ns

from .common_utils import (
    load_h5ad,
    log_section,
    logger
)


def run_prep_binary_scores(h5ad_path, cluster_header, organ, first_author, year):
    """
    Compute binary scores per cluster.

    Loads adata_filtered.h5ad, runs ns.pp.prep_medians() then ns.pp.prep_binary_scores()
    in memory (matching DEMO), saves binary scores csv + pkl.
    """
    log_section("NSForest: Prep Binary Scores")

    cluster_header_safe = cluster_header.replace(" ", "_")
    prefix = f"{organ}_{first_author}_{year}_{cluster_header_safe}"

    adata = load_h5ad(h5ad_path, cluster_header)
    adata_prep = adata.copy()

    logger.info("Running ns.pp.prep_medians()...")
    adata_prep = ns.pp.prep_medians(adata_prep, cluster_header)

    logger.info("Running ns.pp.prep_binary_scores()...")
    adata_prep = ns.pp.prep_binary_scores(adata_prep, cluster_header)

    # gene-by-cluster DataFrame — matching DEMO exactly
    df_binary_scores = adata_prep.varm['binary_scores_' + cluster_header]
    logger.info(f"Binary scores shape: {df_binary_scores.shape}")

    df_binary_scores.to_csv(f"{prefix}_binary_scores.csv")
    df_binary_scores.to_pickle(f"{prefix}_binary_scores.pkl")
    logger.info(f"Saved: {prefix}_binary_scores.csv")
    logger.info(f"Saved: {prefix}_binary_scores.pkl")

    logger.info("Prep binary scores complete!")
