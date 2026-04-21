"""
Run NSForest algorithm to identify marker genes (parallelized by cluster batch).

Corresponds to DEMO_NS-Forest_workflow.py: Section 3 run NSForest()

Loads adata_filtered.h5ad, reads medians and binary_scores CSVs into varm,
then calls nsforesting.NSForest() with cluster_list for parallelization.

Saves:
  {organ}_{first_author}_{journal}_{year}_{cluster_header}_{embedding}_{vid}_results.csv
"""

import pandas as pd
from nsforest import nsforesting

from .common_utils import (
    get_output_prefix,
    load_h5ad,
    log_section,
    logger
)


def run_nsforest(h5ad_path, medians_csv, binary_scores_csv, cluster_header,
                 organ, first_author, journal, year, embedding, dataset_version_id,
                 cluster_list=None, n_trees=1000, n_genes_eval=6):
    """
    Run NSForest for a batch of clusters.

    Parameters
    ----------
    h5ad_path         : Path to adata_filtered.h5ad
    medians_csv       : Path to {prefix}_medians.csv (gene-by-cluster)
    binary_scores_csv : Path to {prefix}_binary_scores.csv (gene-by-cluster)
    cluster_header    : Column name for cell type clusters
    cluster_list      : List of cluster names to process (for parallelization)
    """
    log_section("NSForest: Run NSForest")

    prefix = get_output_prefix( organ, first_author, journal, year, cluster_header, embedding, dataset_version_id )

    # Load filtered adata
    adata_prep = load_h5ad(h5ad_path, cluster_header)
    adata_prep = adata_prep.copy()

    # Load medians and binary scores CSVs (gene-by-cluster, matching DEMO)
    logger.info(f"Loading medians: {medians_csv}")
    df_medians = pd.read_csv(medians_csv, index_col=0)
    logger.info(f"Medians shape: {df_medians.shape}")

    logger.info(f"Loading binary scores: {binary_scores_csv}")
    df_binary_scores = pd.read_csv(binary_scores_csv, index_col=0)
    logger.info(f"Binary scores shape: {df_binary_scores.shape}")

    # Subset adata_prep to positive genes (index of medians CSV)
    adata_prep = adata_prep[:, df_medians.index].copy()

    # Attach to varm — gene-by-cluster, matching DEMO
    adata_prep.varm['medians_' + cluster_header] = df_medians
    adata_prep.varm['binary_scores_' + cluster_header] = df_binary_scores

    # Run NSForest
    if cluster_list:
        logger.info(f"Running NSForest for cluster(s): {cluster_list}")
    else:
        logger.info("Running NSForest for all clusters")

    results = nsforesting.NSForest(
        adata_prep,
        cluster_header,
        cluster_list=cluster_list if cluster_list else [],
        n_trees=n_trees,
        n_genes_eval=n_genes_eval,
        save=False,
        save_supplementary=False,
    )

    logger.info(f"NSForest results shape: {results.shape}")

    # Save partial results — unique filename per batch
    if cluster_list:
        cluster_safe = cluster_list[0].replace(' ', '_').replace('/', '-')
        output_csv = f"{prefix}_results_{cluster_safe}.csv"
    else:
        output_csv = f"{prefix}_results.csv"

    results.to_csv(output_csv, index=False)
    logger.info(f"Saved: {output_csv}")
    
    logger.info("NSForest complete!")
