"""
Compute median expression and binary scores per cluster.

Runs ns.pp.prep_medians() and ns.pp.prep_binary_scores() sequentially on the
same adata_prep object, matching DEMO_NS-forest_workflow.py Section 3 exactly:

    adata_prep = adata.copy()
    adata_prep = ns.pp.prep_medians(adata_prep, cluster_header)
    adata_prep = ns.pp.prep_binary_scores(adata_prep, cluster_header)

varm keys are renamed to h5py-safe names before writing adata_prep.h5ad
to avoid TypeError when cluster names contain '/' (e.g. 'VSMC/P').
run_nsforest.py renames them back before calling nsforesting.NSForest().
"""

import nsforest as ns

from .common_utils import (
    create_output_dir,
    load_h5ad,
    log_section,
    logger
)

# Sentinel used to make varm keys h5py-safe.
# Must match the sentinel used in run_nsforest.py.
_SLASH = "__SLASH__"


def _safe_key(key):
    """Replace '/' with _SLASH so h5py can write the key."""
    return key.replace("/", _SLASH)


def run_prep_medians_binary_scores(h5ad_path, cluster_header, organ, first_author, year):
    """
    Compute medians and binary scores, then save adata_prep.h5ad + CSVs.

    Saves:
    - adata_prep.h5ad          (varm keys sanitized for h5py)
    - {cluster_header}_medians.csv
    - {cluster_header}_binary_scores.csv
    """
    log_section("NSForest: Prep Medians + Binary Scores")

    output_folder = create_output_dir(organ, first_author, year)
    outputfilename_prefix = cluster_header.replace(" ", "_")

    # Load and copy — matching DEMO exactly
    adata = load_h5ad(h5ad_path, cluster_header)
    adata_prep = adata.copy()

    # Step 1: medians
    logger.info("Running ns.pp.prep_medians()...")
    adata_prep = ns.pp.prep_medians(adata_prep, cluster_header)

    # Step 2: binary scores (requires medians_ already in varm)
    logger.info("Running ns.pp.prep_binary_scores()...")
    adata_prep = ns.pp.prep_binary_scores(adata_prep, cluster_header)

    medians_key       = 'medians_' + cluster_header
    binary_scores_key = 'binary_scores_' + cluster_header

    # Save CSVs
    df_medians       = adata_prep.varm[medians_key].T
    df_binary_scores = adata_prep.varm[binary_scores_key].T

    logger.info(f"Median matrix shape: {df_medians.shape}")
    medians_csv = f"{output_folder}/{outputfilename_prefix}_medians.csv"
    df_medians.to_csv(medians_csv)
    logger.info(f"Saved: {medians_csv}")

    logger.info(f"Binary scores shape: {df_binary_scores.shape}")
    binary_scores_csv = f"{output_folder}/{outputfilename_prefix}_binary_scores.csv"
    df_binary_scores.to_csv(binary_scores_csv)
    logger.info(f"Saved: {binary_scores_csv}")

    # Rename varm keys to h5py-safe names before writing h5ad.
    # '/' in cluster names (e.g. 'VSMC/P') causes h5py TypeError.
    # run_nsforest.py reverses this rename before calling nsforesting.NSForest().
    safe_medians_key       = _safe_key(medians_key)
    safe_binary_scores_key = _safe_key(binary_scores_key)

    adata_prep.varm[safe_medians_key]       = adata_prep.varm[medians_key]
    adata_prep.varm[safe_binary_scores_key] = adata_prep.varm[binary_scores_key]
    del adata_prep.varm[medians_key]
    del adata_prep.varm[binary_scores_key]

    adata_prep_path = f"{output_folder}/adata_prep.h5ad"
    adata_prep.write_h5ad(adata_prep_path)
    logger.info(f"Saved: adata_prep.h5ad")

    logger.info("Prep medians + binary scores complete!")
