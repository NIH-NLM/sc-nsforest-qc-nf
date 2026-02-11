"""
Filter adata by disease status, tissue, and minimum cluster size.

Creates before/after dendrograms and cluster statistics to show filtering effects.
"""

import pandas as pd
import nsforest as ns

from .common_utils import (
    create_output_dir,
    load_h5ad,
    log_section,
    logger
)


def create_stats_before_filter(adata, cluster_header, output_folder, outputfilename_prefix):
    """Create dendrogram and stats BEFORE filtering."""
    logger.info("Creating BEFORE FILTER statistics...")
    
    n_clusters = adata.obs[cluster_header].nunique()
    logger.info(f"Before filter - Total cells: {adata.n_obs}, Clusters: {n_clusters}")
    
    fig_width = int(n_clusters / 5)
    fig_height = max([2, int(max([len(z) for z in adata.obs[cluster_header].unique()]) / 30) + 1])
    
    ns.pp.dendrogram(
        adata, cluster_header, figsize=(fig_width, fig_height),
        tl_kwargs={'optimal_ordering': True}, save="svg",
        output_folder=output_folder,
        outputfilename_suffix=outputfilename_prefix + "_before_filter"
    )
    
    # Save cluster sizes with proper column names
    cluster_counts = adata.obs[cluster_header].value_counts()
    df_cluster_sizes = pd.DataFrame({
        'cluster': cluster_counts.index.tolist(),
        'count': cluster_counts.values
    })
    df_cluster_sizes.to_csv(output_folder + "/" + outputfilename_prefix + "_cluster_sizes_before_filter.csv", index=False)
    
    cluster_order = [x.strip() for x in adata.uns["dendrogram_" + cluster_header]['categories_ordered']]
    pd.DataFrame({'cluster_order': cluster_order}).to_csv(
        output_folder + "/" + outputfilename_prefix + "_cluster_order_before_filter.csv", 
        index=False
    )
    
    df_summary = pd.DataFrame({
        'n_obs': [adata.n_obs], 'n_vars': [adata.n_vars], 'n_clusters': [n_clusters]
    })
    df_summary.to_csv(output_folder + "/" + outputfilename_prefix + "_summary_before_filter.csv", index=False)
    
    logger.info("Before filter statistics saved")


def filter_by_disease(adata, disease_column='disease', filter_normal=False):
    """Filter cells by disease status."""
    if not filter_normal:
        logger.info("Filter normal = False, keeping all disease states")
        return adata
    
    if disease_column not in adata.obs.columns:
        logger.warning(f"Disease column '{disease_column}' not found in adata.obs")
        logger.warning("Skipping disease filtering")
        return adata
    
    logger.info(f"Filtering to normal cells only (disease column: {disease_column})")
    
    mask = adata.obs[disease_column].str.lower().str.contains('normal', na=False)
    n_before = adata.n_obs
    adata_filtered = adata[mask].copy()
    n_after = adata_filtered.n_obs
    
    logger.info(f"Disease filtering: {n_before} → {n_after} cells ({n_before - n_after} removed)")
    
    return adata_filtered


def filter_by_tissue(adata, tissue_column='tissue', target_tissue=None):
    """Filter cells by tissue."""
    if target_tissue is None:
        logger.info("No tissue filter specified, keeping all tissues")
        return adata
    
    if tissue_column not in adata.obs.columns:
        logger.warning(f"Tissue column '{tissue_column}' not found in adata.obs")
        logger.warning("Skipping tissue filtering")
        return adata
    
    logger.info(f"Filtering to tissue: {target_tissue} (tissue column: {tissue_column})")
    
    mask = adata.obs[tissue_column].str.lower().str.contains(target_tissue.lower(), na=False)
    n_before = adata.n_obs
    adata_filtered = adata[mask].copy()
    n_after = adata_filtered.n_obs
    
    logger.info(f"Tissue filtering: {n_before} → {n_after} cells ({n_before - n_after} removed)")
    
    return adata_filtered


def filter_by_min_cluster_size(adata, cluster_header, min_size=5):
    """Remove clusters with fewer than min_size cells."""
    logger.info(f"Filtering clusters with < {min_size} cells")
    
    cluster_counts = adata.obs[cluster_header].value_counts()
    clusters_to_keep = cluster_counts[cluster_counts >= min_size].index
    
    n_clusters_before = adata.obs[cluster_header].nunique()
    n_clusters_after = len(clusters_to_keep)
    removed_clusters = n_clusters_before - n_clusters_after
    
    if removed_clusters > 0:
        logger.info(f"Removing {removed_clusters} clusters with < {min_size} cells:")
        small_clusters = cluster_counts[cluster_counts < min_size]
        for cluster, count in small_clusters.items():
            logger.info(f"  - {cluster}: {count} cells")
    
    mask = adata.obs[cluster_header].isin(clusters_to_keep)
    n_cells_before = adata.n_obs
    adata_filtered = adata[mask].copy()
    n_cells_after = adata_filtered.n_obs
    
    logger.info(f"Min cluster size filtering: {n_cells_before} cells, {n_clusters_before} clusters "
                f"→ {n_cells_after} cells, {n_clusters_after} clusters")
    
    return adata_filtered


def run_filter_adata(h5ad_path, cluster_header, organ, first_author, year,
                     filter_normal=False, tissue=None, disease_column='disease', 
                     tissue_column='tissue', min_cluster_size=5):
    """
    Filter adata and create before/after statistics.
    
    Filtering steps:
    1. Filter by disease (if filter_normal=True, keep only 'normal' cells)
    2. Filter by tissue (if specified)
    3. Remove clusters with < min_cluster_size cells
    """
    log_section("NSForest: Filter AnnData")
    
    output_folder = create_output_dir(organ, first_author, year) + "/"
    outputfilename_prefix = cluster_header
    
    logger.info(f"Loading: {h5ad_path}")
    adata = load_h5ad(h5ad_path, cluster_header)
    logger.info(f"Original data: {adata.n_obs} cells, {adata.n_vars} genes, "
                f"{adata.obs[cluster_header].nunique()} clusters")
    
    create_stats_before_filter(adata, cluster_header, output_folder, outputfilename_prefix)
    
    logger.info("\n=== Applying Filters ===")
    
    adata_filtered = filter_by_disease(adata, disease_column, filter_normal)
    adata_filtered = filter_by_tissue(adata_filtered, tissue_column, tissue)
    adata_filtered = filter_by_min_cluster_size(adata_filtered, cluster_header, min_cluster_size)
    
    logger.info(f"\nFinal filtered data: {adata_filtered.n_obs} cells, {adata_filtered.n_vars} genes, "
                f"{adata_filtered.obs[cluster_header].nunique()} clusters")
    
    logger.info("\n=== Creating AFTER FILTER statistics ===")
    
    n_clusters = adata_filtered.obs[cluster_header].nunique()
    fig_width = int(n_clusters / 5)
    fig_height = max([2, int(max([len(z) for z in adata_filtered.obs[cluster_header].unique()]) / 30) + 1])
    
    ns.pp.dendrogram(
        adata_filtered, cluster_header, figsize=(fig_width, fig_height),
        tl_kwargs={'optimal_ordering': True}, save="svg",
        output_folder=output_folder, outputfilename_suffix=outputfilename_prefix
    )
    
    # Save cluster sizes with proper column names
    cluster_counts = adata_filtered.obs[cluster_header].value_counts()
    df_cluster_sizes = pd.DataFrame({
        'cluster': cluster_counts.index.tolist(),
        'count': cluster_counts.values
    })
    df_cluster_sizes.to_csv(output_folder + outputfilename_prefix + "_cluster_sizes.csv", index=False)
    
    cluster_order = [x.strip() for x in adata_filtered.uns["dendrogram_" + cluster_header]['categories_ordered']]
    pd.DataFrame({'cluster_order': cluster_order}).to_csv(
        output_folder + outputfilename_prefix + "_cluster_order.csv", 
        index=False
    )
    
    df_normal = pd.DataFrame({
        'n_obs': [adata_filtered.n_obs], 'n_vars': [adata_filtered.n_vars], 'n_clusters': [n_clusters]
    })
    df_normal.to_csv(output_folder + outputfilename_prefix + "_summary_normal.csv", index=False)
    
    filtered_h5ad_path = output_folder + "adata_filtered.h5ad"
    adata_filtered.write_h5ad(filtered_h5ad_path)
    logger.info(f"Saved filtered h5ad: {filtered_h5ad_path}")
    
    logger.info("\nFilter complete!")
