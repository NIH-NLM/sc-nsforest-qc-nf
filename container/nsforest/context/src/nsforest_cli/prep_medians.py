"""
Compute median expression per cluster.

Corresponds to DEMO_NS-Forest_workflow.py: Section 3 prep
Uses ns.pp.prep_medians() to filter positive genes and compute medians.

Saves:
  {organ}_{first_author}_{journal}_{year}_{cluster_header}_{embedding}_{vid}_medians.csv
  {organ}_{first_author}_{journal}_{year}_{cluster_header}_{embedding}_{vid}_medians.pkl
"""

import nsforest as ns

from .common_utils import (
    get_output_prefix,
    load_h5ad,
    log_section,
    logger
)


def run_prep_medians(h5ad_path, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id):
    """
    Compute median expression per cluster.

    Loads adata_filtered.h5ad, runs ns.pp.prep_medians(), saves medians csv + pkl.
    """
    log_section("NSForest: Prep Medians")

    prefix = get_output_prefix( organ, first_author, journal, year, cluster_header, embedding, dataset_version_id )

    adata = load_h5ad(h5ad_path, cluster_header)
    adata_prep = adata.copy()

    logger.info("Running ns.pp.prep_medians()...")
    adata_prep = ns.pp.prep_medians(adata_prep, cluster_header)

    # gene-by-cluster DataFrame — matching DEMO exactly
    df_medians = adata_prep.varm['medians_' + cluster_header]
    logger.info(f"Medians shape: {df_medians.shape}")

    df_medians.to_csv(f"{prefix}_medians.csv")
    df_medians.to_pickle(f"{prefix}_medians.pkl")
    logger.info(f"Saved: {prefix}_medians.csv")
    logger.info(f"Saved: {prefix}_medians.pkl")

    logger.info("Prep medians complete!")
