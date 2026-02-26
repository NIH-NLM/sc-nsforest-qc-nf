"""
Filter adata by disease, tissue, age, and minimum cluster size.

Creates before/after dendrograms and cluster statistics to show filtering effects.

All filtering uses ontology term ID columns only - no text matching:
  - Tissue  : tissue_ontology_term_id   IN  UBERON obo_ids  (uberon_json)
  - Disease : disease_ontology_term_id  IN  PATO/MONDO obo_ids  (disease_json)
  - Age     : development_stage_ontology_term_id resolved via HsapDv ages (hsapdv_json)
"""

import json
import pandas as pd
import nsforest as ns

from .common_utils import (
    create_output_dir,
    load_h5ad,
    log_section,
    logger
)


# =============================================================================
# Ontology JSON loaders
# =============================================================================

def load_obo_ids(json_path: str, label: str) -> set:
    """Load obo_ids set from a resolve_uberon / resolve_disease JSON file."""
    with open(json_path) as f:
        data = json.load(f)
    obo_ids = set(data["obo_ids"])
    roots   = [t["label"] for t in data["root_terms"]]
    logger.info(f"  Loaded {label} JSON : {json_path}")
    logger.info(f"  Root terms          : {', '.join(roots)}")
    logger.info(f"  Total obo_ids       : {len(obo_ids):,}")
    return obo_ids


# =============================================================================
# Before-filter statistics
# =============================================================================

def create_stats_before_filter(adata, cluster_header, output_folder, outputfilename_prefix):
    """Create dendrogram and stats BEFORE filtering."""
    logger.info("Creating BEFORE FILTER statistics...")

    n_clusters = adata.obs[cluster_header].nunique()
    logger.info(f"Before filter - Total cells: {adata.n_obs}, Clusters: {n_clusters}")

    fig_width  = int(n_clusters / 5)
    fig_height = max([2, int(max([len(z) for z in adata.obs[cluster_header].unique()]) / 30) + 1])

    ns.pp.dendrogram(
        adata, cluster_header, figsize=(fig_width, fig_height),
        tl_kwargs={"optimal_ordering": True}, save="svg",
        output_folder=output_folder,
        outputfilename_suffix=outputfilename_prefix + "_before_filter"
    )

    cluster_counts   = adata.obs[cluster_header].value_counts()
    df_cluster_sizes = pd.DataFrame({
        "cluster": cluster_counts.index.tolist(),
        "count":   cluster_counts.values
    })
    df_cluster_sizes.to_csv(
        output_folder + "/" + outputfilename_prefix + "_cluster_sizes_before_filter.csv",
        index=False
    )

    cluster_order = [
        x.strip() for x in adata.uns["dendrogram_" + cluster_header]["categories_ordered"]
    ]
    pd.DataFrame({"cluster_order": cluster_order}).to_csv(
        output_folder + "/" + outputfilename_prefix + "_cluster_order_before_filter.csv",
        index=False
    )

    pd.DataFrame({
        "n_obs": [adata.n_obs], "n_vars": [adata.n_vars], "n_clusters": [n_clusters]
    }).to_csv(
        output_folder + "/" + outputfilename_prefix + "_summary_before_filter.csv",
        index=False
    )

    logger.info("Before filter statistics saved")


# =============================================================================
# Filter steps - ontology ID columns only
# =============================================================================

def filter_by_tissue(adata, uberon_json=None):
    """Filter cells using tissue_ontology_term_id matched against UBERON obo_ids.

    Args:
        adata:       AnnData object
        uberon_json: Path to UBERON JSON from cellxgene-harvester resolve-uberon.
                     Skipped if None.
    """
    if not uberon_json:
        logger.info("No --uberon file provided - skipping tissue filter")
        return adata

    id_col = "tissue_ontology_term_id"
    if id_col not in adata.obs.columns:
        raise ValueError(
            f"'{id_col}' not found in adata.obs. "
            "This column is required for UBERON-based tissue filtering. "
            "Ensure the h5ad file originates from CellxGene."
        )

    obo_ids  = load_obo_ids(uberon_json, "UBERON")
    mask     = adata.obs[id_col].isin(obo_ids)
    n_before = adata.n_obs
    adata    = adata[mask].copy()

    logger.info(f"Tissue filter (UBERON IDs): {n_before} -> {adata.n_obs} cells "
                f"({n_before - adata.n_obs} removed)")
    return adata


def filter_by_disease(adata, disease_json=None, filter_normal=False):
    """Filter cells using disease_ontology_term_id matched against disease obo_ids.

    Args:
        adata:         AnnData object
        disease_json:  Path to disease JSON from cellxgene-harvester resolve-disease.
                       Skipped if None or filter_normal=False.
        filter_normal: If False, no disease filtering is applied.
    """
    if not filter_normal:
        logger.info("filter_normal=False - keeping all disease states")
        return adata

    if not disease_json:
        logger.warning("filter_normal=True but no --disease file provided - skipping disease filter")
        return adata

    id_col = "disease_ontology_term_id"
    if id_col not in adata.obs.columns:
        raise ValueError(
            f"'{id_col}' not found in adata.obs. "
            "This column is required for disease filtering. "
            "Ensure the h5ad file originates from CellxGene."
        )

    obo_ids  = load_obo_ids(disease_json, "disease")
    mask     = adata.obs[id_col].isin(obo_ids)
    n_before = adata.n_obs
    adata    = adata[mask].copy()

    logger.info(f"Disease filter (ontology IDs): {n_before} -> {adata.n_obs} cells "
                f"({n_before - adata.n_obs} removed)")
    return adata


def filter_by_age(adata, hsapdv_json=None, filter_normal=False):
    """Filter cells using development_stage_ontology_term_id matched against HsapDv obo_ids.

    The age threshold is encoded in the JSON at resolve time
    (cellxgene-harvester resolve-hsapdv --min-age N).
    No numeric comparison here — identical pattern to filter_by_tissue and filter_by_disease.

    Args:
        adata:         AnnData object
        hsapdv_json:   Path to HsapDv JSON from cellxgene-harvester resolve-hsapdv --min-age N.
                       Skipped if None or filter_normal=False.
        filter_normal: If False, no age filtering is applied.
    """
    if not filter_normal:
        logger.info("filter_normal=False - skipping age filter")
        return adata

    if not hsapdv_json:
        logger.warning("filter_normal=True but no --hsapdv file provided - skipping age filter")
        return adata

    id_col = "development_stage_ontology_term_id"
    if id_col not in adata.obs.columns:
        raise ValueError(
            f"'{id_col}' not found in adata.obs. "
            "This column is required for age filtering. "
            "Ensure the h5ad file originates from CellxGene."
        )

    obo_ids  = load_obo_ids(hsapdv_json, "HsapDv")
    mask     = adata.obs[id_col].isin(obo_ids)
    n_before = adata.n_obs

    # Log every term ID present in the data with its kept/excluded status and cell count
    all_counts = adata.obs[id_col].value_counts()
    kept_terms     = [(tid, cnt) for tid, cnt in all_counts.items() if tid in obo_ids]
    excluded_terms = [(tid, cnt) for tid, cnt in all_counts.items() if tid not in obo_ids]

    logger.info(f"  Kept     ({len(kept_terms)} terms, "
                f"{sum(c for _,c in kept_terms):,} cells):")
    for term_id, count in kept_terms:
        logger.info(f"    KEPT     {term_id}: {count:,} cells")

    if excluded_terms:
        logger.info(f"  Excluded ({len(excluded_terms)} terms, "
                    f"{sum(c for _,c in excluded_terms):,} cells):")
        for term_id, count in excluded_terms:
            logger.info(f"    EXCLUDED {term_id}: {count:,} cells")
    else:
        logger.info("  No cells excluded by age filter")

    adata = adata[mask].copy()
    logger.info(f"Age filter (HsapDv IDs): {n_before} -> {adata.n_obs} cells "
                f"({n_before - adata.n_obs} removed)")
    return adata


def filter_by_min_cluster_size(adata, cluster_header, min_size=5):
    """Remove clusters with fewer than min_size cells."""
    logger.info(f"Filtering clusters with < {min_size} cells")

    cluster_counts    = adata.obs[cluster_header].value_counts()
    clusters_to_keep  = cluster_counts[cluster_counts >= min_size].index
    n_clusters_before = adata.obs[cluster_header].nunique()
    n_clusters_after  = len(clusters_to_keep)
    removed_clusters  = n_clusters_before - n_clusters_after

    if removed_clusters > 0:
        logger.info(f"Removing {removed_clusters} clusters with < {min_size} cells:")
        for cluster, count in cluster_counts[cluster_counts < min_size].items():
            logger.info(f"  - {cluster}: {count} cells")

    n_cells_before = adata.n_obs
    adata          = adata[adata.obs[cluster_header].isin(clusters_to_keep)].copy()

    logger.info(f"Min cluster size filter: {n_cells_before} cells, {n_clusters_before} clusters "
                f"-> {adata.n_obs} cells, {n_clusters_after} clusters")
    return adata


# =============================================================================
# Main entry point
# =============================================================================

def run_filter_adata(h5ad_path, cluster_header, organ, first_author, year,
                     filter_normal=False,
                     uberon_json=None,
                     disease_json=None,
                     hsapdv_json=None,
                     min_cluster_size=5):
    """
    Filter adata and create before/after statistics.

    Filtering steps (all via ontology term ID columns - no text matching):
      1. Tissue  - tissue_ontology_term_id   IN UBERON obo_ids  (uberon_json)
      2. Disease - disease_ontology_term_id  IN disease obo_ids (disease_json)
      3. Age     - development_stage_ontology_term_id IN HsapDv obo_ids (hsapdv_json)
                   Age threshold is encoded in hsapdv_json at resolve time:
                   cellxgene-harvester resolve-hsapdv --min-age 15
      4. Min cluster size

    Args:
        h5ad_path:        Path to input h5ad file
        cluster_header:   Column name for cell type clusters
        organ:            Organ/tissue label (used for output directory)
        first_author:     First author (used for output directory)
        year:             Publication year (used for output directory)
        filter_normal:    If True, apply tissue + disease + age filters
        uberon_json:      Path to UBERON JSON from resolve-uberon
        disease_json:     Path to disease JSON from resolve-disease
        hsapdv_json:      Path to HsapDv JSON from resolve-hsapdv --min-age N
        min_cluster_size: Minimum cells per cluster (default 5)
    """
    log_section("NSForest: Filter AnnData")

    output_folder         = create_output_dir(organ, first_author, year) + "/"
    outputfilename_prefix = cluster_header

    logger.info(f"Loading: {h5ad_path}")
    adata = load_h5ad(h5ad_path, cluster_header)
    logger.info(f"Original data: {adata.n_obs} cells, {adata.n_vars} genes, "
                f"{adata.obs[cluster_header].nunique()} clusters")

    create_stats_before_filter(adata, cluster_header, output_folder, outputfilename_prefix)

    logger.info("\n=== Applying Filters ===")

    if filter_normal:
        logger.info("\n[1/4] Tissue filter")
        adata = filter_by_tissue(adata, uberon_json)

        logger.info("\n[2/4] Disease filter")
        adata = filter_by_disease(adata, disease_json, filter_normal)

        logger.info("\n[3/4] Age filter")
        adata = filter_by_age(adata, hsapdv_json, filter_normal)
    else:
        logger.info("filter_normal=False - skipping tissue, disease, and age filters")

    logger.info("\n[4/4] Min cluster size filter")
    adata = filter_by_min_cluster_size(adata, cluster_header, min_cluster_size)

    logger.info(f"\nFinal filtered data: {adata.n_obs} cells, {adata.n_vars} genes, "
                f"{adata.obs[cluster_header].nunique()} clusters")

    logger.info("\n=== Creating AFTER FILTER statistics ===")

    n_clusters = adata.obs[cluster_header].nunique()
    fig_width  = int(n_clusters / 5)
    fig_height = max([2, int(max([len(z) for z in adata.obs[cluster_header].unique()]) / 30) + 1])

    ns.pp.dendrogram(
        adata, cluster_header, figsize=(fig_width, fig_height),
        tl_kwargs={"optimal_ordering": True}, save="svg",
        output_folder=output_folder, outputfilename_suffix=outputfilename_prefix
    )

    cluster_counts   = adata.obs[cluster_header].value_counts()
    df_cluster_sizes = pd.DataFrame({
        "cluster": cluster_counts.index.tolist(),
        "count":   cluster_counts.values
    })
    df_cluster_sizes.to_csv(
        output_folder + outputfilename_prefix + "_cluster_sizes.csv", index=False
    )

    cluster_order = [
        x.strip() for x in adata.uns["dendrogram_" + cluster_header]["categories_ordered"]
    ]
    pd.DataFrame({"cluster_order": cluster_order}).to_csv(
        output_folder + outputfilename_prefix + "_cluster_order.csv", index=False
    )

    pd.DataFrame({
        "n_obs": [adata.n_obs], "n_vars": [adata.n_vars], "n_clusters": [n_clusters]
    }).to_csv(
        output_folder + outputfilename_prefix + "_summary_normal.csv", index=False
    )

    filtered_h5ad_path = output_folder + "adata_filtered.h5ad"
    adata.write_h5ad(filtered_h5ad_path)
    logger.info(f"Saved filtered h5ad: {filtered_h5ad_path}")

    logger.info("\nFilter complete!")
