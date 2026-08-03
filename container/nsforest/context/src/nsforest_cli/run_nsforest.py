"""
Run NSForest algorithm to identify marker genes (parallelized by cluster batch).

Corresponds to DEMO_NS-Forest_workflow.py: Section 3 run NSForest()

Loads adata_filtered.h5ad, reads medians and binary_scores CSVs into varm,
then calls nsforesting.NSForest() with cluster_list for parallelization.

Saves:
  results_{organ}_{first_author}_{journal}_{year}_{cluster_header}_{embedding}_{vid}.csv
"""
import csv
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
                 cluster_list=None, n_trees=1000, n_genes_eval=6,
                 max_cells_per_cluster=0, seed=42):
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

    # Load filtered adata (no .copy(): the positive-gene subset below makes a fresh object,
    # so an extra full-matrix copy here just doubled peak RAM on large datasets).
    adata_prep = load_h5ad(h5ad_path, cluster_header)

    # Load medians and binary scores CSVs (gene-by-cluster, matching DEMO)
    logger.info(f"Loading medians: {medians_csv}")
    df_medians = pd.read_csv(medians_csv, index_col=0)
    logger.info(f"Medians shape: {df_medians.shape}")

    logger.info(f"Loading binary scores: {binary_scores_csv}")
    df_binary_scores = pd.read_csv(binary_scores_csv, index_col=0)
    logger.info(f"Binary scores shape: {df_binary_scores.shape}")

    # Subset adata_prep to positive genes (index of medians CSV)
    adata_prep = adata_prep[:, df_medians.index].copy()

    # Diagnose any mismatch between df_medians.index and adata_prep.var_names
    same_len = len(df_medians.index) == len(adata_prep.var_names)
    if not same_len:
        missing = sorted(set(df_medians.index) - set(adata_prep.var_names))
        logger.warning(
            f"Gene count mismatch — df_medians: {len(df_medians.index)}, "
            f"adata.var_names: {len(adata_prep.var_names)}. "
            f"Missing from adata ({len(missing)}): {missing[:5]}..."
        )
    elif not (df_medians.index == adata_prep.var_names).all():
        logger.warning(
            f"Order/dtype mismatch — df_medians dups: {df_medians.index.has_duplicates}, "
            f"var_names dups: {adata_prep.var_names.has_duplicates}, "
            f"df_medians dtype: {df_medians.index.dtype}, "
            f"var_names dtype: {adata_prep.var_names.dtype}"
        )

    # Defensive — force exact alignment to adata_prep.var_names so varm validation passes.
    # If reindex introduces NaN rows, the underlying mismatch is real (will fail loudly below).
    df_medians       = df_medians.reindex(adata_prep.var_names)
    df_binary_scores = df_binary_scores.reindex(adata_prep.var_names)

    if df_medians.isna().any().any() or df_binary_scores.isna().any().any():
        n_bad = df_medians.isna().any(axis=1).sum()
        raise ValueError(
            f"reindex left {n_bad} rows as NaN — df_medians.index doesn't fully cover "
            f"adata_prep.var_names. Investigate var_names duplicates or stale CSV."
        )

    # Attach to varm — gene-by-cluster, matching DEMO
    adata_prep.varm['medians_' + cluster_header] = df_medians
    adata_prep.varm['binary_scores_' + cluster_header] = df_binary_scores

    # Optional: cap cells per cluster for the RandomForest / evaluation (memory + time at scale).
    # Applied AFTER attaching full-data medians/binary_scores to varm, so candidate-gene
    # selection stays full-data; only the RF training + F-beta evaluation see the subsample.
    # Stratified + seeded → reproducible; clusters smaller than the cap keep all their cells.
    if max_cells_per_cluster and max_cells_per_cluster > 0:
        import numpy as np
        rng = np.random.default_rng(seed)
        groups = adata_prep.obs.groupby(cluster_header, observed=True).indices
        keep = []
        for cl, idx in groups.items():
            idx = np.asarray(idx)
            if len(idx) > max_cells_per_cluster:
                idx = rng.choice(idx, size=max_cells_per_cluster, replace=False)
            keep.append(idx)
        keep_idx = np.sort(np.concatenate(keep))
        n_before = adata_prep.n_obs
        adata_prep = adata_prep[keep_idx].copy()
        logger.info(
            f"Subsampled cells for NSForest: {n_before} -> {adata_prep.n_obs} "
            f"(max {max_cells_per_cluster}/cluster, seed={seed}). "
            f"Medians/binary scores remain full-data."
        )

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
        # in case there is a problem with a stray quote - clean it up before output
        cluster_safe = cluster_list[0].replace('"', '').replace("'", '').replace(' ', '_').replace('/', '-')
        output_csv = f"results_{cluster_safe}_{prefix}.csv"
    else:
        output_csv = f"results_{prefix}.csv"

    results.to_csv(output_csv, index=False)
    logger.info(f"Saved: {output_csv}")
    
    logger.info("NSForest complete!")
