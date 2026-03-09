"""
Compute median expression and binary scores per cluster.

Runs ns.pp.prep_medians() and ns.pp.prep_binary_scores() sequentially on the
same adata_prep object, matching DEMO_NS-forest_workflow.py Section 3 exactly:

    adata_prep = adata.copy()
    adata_prep = ns.pp.prep_medians(adata_prep, cluster_header)
    adata_prep = ns.pp.prep_binary_scores(adata_prep, cluster_header)

Both varm entries are extracted and deleted before writing adata_prep.h5ad to
avoid h5py TypeError when cluster names contain '/' (e.g. 'VSMC/P').
"""

import nsforest as ns

from .common_utils import (
    create_output_dir,
    load_h5ad,
    log_section,
    logger
)


def run_prep_medians_binary_scores(h5ad_path, cluster_header, organ, first_author, year):
    """
    Compute medians and binary scores, then save adata_prep.h5ad + CSVs.

    Saves:
    - adata_prep.h5ad          (positive genes only, varm cleared of medians/binary scores)
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

    # Extract both dataframes from varm before writing h5ad.
    # Cluster names (e.g. 'VSMC/P') contain '/' which h5py treats as a path
    # separator and raises TypeError when writing as dataset keys.
    df_medians = adata_prep.varm['medians_' + cluster_header].T
    del adata_prep.varm['medians_' + cluster_header]

    df_binary_scores = adata_prep.varm['binary_scores_' + cluster_header].T
    del adata_prep.varm['binary_scores_' + cluster_header]

    # Write h5ad — now safe, no cluster names in varm keys
    adata_prep_path = f"{output_folder}/adata_prep.h5ad"
    adata_prep.write_h5ad(adata_prep_path)
    logger.info(f"Saved: adata_prep.h5ad")

    # Save medians CSV
    logger.info(f"Median matrix shape: {df_medians.shape}")
    medians_csv = f"{output_folder}/{outputfilename_prefix}_medians.csv"
    df_medians.to_csv(medians_csv)
    logger.info(f"Saved: {medians_csv}")

    # Save binary scores CSV
    logger.info(f"Binary scores shape: {df_binary_scores.shape}")
    binary_scores_csv = f"{output_folder}/{outputfilename_prefix}_binary_scores.csv"
    df_binary_scores.to_csv(binary_scores_csv)
    logger.info(f"Saved: {binary_scores_csv}")

    logger.info("Prep medians + binary scores complete!")
