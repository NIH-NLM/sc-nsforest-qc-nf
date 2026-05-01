"""
Filter adata by disease, tissue, age, and minimum cluster size.

Creates before/after dendrograms and cluster statistics to show filtering effects.

All filtering uses ontology term ID columns only - no text matching:
  - Tissue  : tissue_ontology_term_id   IN  UBERON obo_ids  (uberon_json)
  - Disease : disease_ontology_term_id  IN  PATO/MONDO obo_ids  (disease_json)
  - Age     : development_stage_ontology_term_id IN HsapDv obo_ids (hsapdv_json)
              Age threshold is encoded in the JSON at resolve time.

Additionally, add gene_symbols to all the ensgs and store the new variable in the adata object
"""

# Force non-interactive backend before any other imports that might trigger a display.
# Required for headless execution in Nextflow / Docker containers.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import json
import pandas as pd
import scanpy as sc
import nsforest as ns

from .common_utils import (
    get_output_prefix,
    load_h5ad,
    log_section,
    logger
)
from .gene_mapping_utils import (
    load_gene_mapping,
    add_gene_symbols_to_adata
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



def load_obo_labels(json_path: str) -> dict:
    """Load obo_id -> label mapping from any resolve JSON file.

    Used to annotate log output with human-readable labels alongside term IDs.
    """
    with open(json_path) as f:
        data = json.load(f)
    return {t["obo_id"]: t["label"] for t in data["terms"] if t.get("obo_id")}

# =============================================================================
# Before-filter statistics
# =============================================================================

def create_stats_before_filter(adata, cluster_header, prefix, embedding="X_pca"):
    """Create dendrogram and stats BEFORE filtering."""
    logger.info("Creating BEFORE FILTER statistics...")

    n_clusters = adata.obs[cluster_header].nunique()
    logger.info(f"Before filter - Total cells: {adata.n_obs}, Clusters: {n_clusters}")

    try:
        # Ensure X_pca exists — ns.pp.dendrogram hardcodes use_rep="X_pca"
        # but its guard is buggy when other keys (e.g. X_umap) are in obsm.
        if "X_pca" not in adata.obsm:
            logger.info("X_pca not found in obsm — computing PCA for dendrogram")
            sc.pp.pca(adata)
        # Render the figure via nsforest/scanpy but save it ourselves.
        # ns.pp.dendrogram(save="svg") routes through scanpy's savefig_or_show,
        # which in scanpy 1.9.6 + matplotlib 3.8.0 crashes with
        # "'bool' object has no attribute 'write'" inside PIL's PNG backend.
        ns.pp.dendrogram(
            adata, cluster_header,
            tl_kwargs={"optimal_ordering": True},
            pl_kwargs={"show": False},
            save=False,
        )
        plt.savefig(
            f"{prefix}_dendrogram_before_filter.svg",
            format="svg", bbox_inches="tight",
        )
        plt.close()

    except ValueError as e:
        if "negative distances" in str(e):
            logger.warning(f"Optimal ordering failed — retrying without: {e}")
            ns.pp.dendrogram(
                adata, cluster_header,
                tl_kwargs={"optimal_ordering": False},
                pl_kwargs={"show": False},
                save=False,
            )
            plt.savefig(f"{prefix}_dendrogram_before_filter.svg", format="svg", bbox_inches="tight")
            plt.close()
        else:
            raise

    cluster_counts   = adata.obs[cluster_header].value_counts()
    df_cluster_sizes = pd.DataFrame({
        "cluster": cluster_counts.index.tolist(),
        "count":   cluster_counts.values
    })
    df_cluster_sizes.to_csv(f"{prefix}_cluster_sizes_before_filter.csv", index=False)

    cluster_order = [
        x.strip() for x in adata.uns["dendrogram_" + cluster_header]["categories_ordered"]
    ]
    pd.DataFrame({"cluster_order": cluster_order}).to_csv(
        f"{prefix}_cluster_order_before_filter.csv", index=False
    )

    pd.DataFrame({
        "n_obs": [adata.n_obs], "n_vars": [adata.n_vars], "n_clusters": [n_clusters]
    }).to_csv(f"{prefix}_summary_before_filter.csv", index=False)

    logger.info("Before filter statistics saved")


# =============================================================================
# Filter steps - ontology ID columns only
# =============================================================================

def filter_by_tissue(adata, uberon_json=None, row_ids=None):
    """Filter cells using tissue_ontology_term_id.

    Uses row_ids (pipe-separated IDs from the CSV row's tissue_ontology_term_id
    column) when provided — these are the pre-computed intersection of dataset
    tissues with the organ UBERON JSON. Falls back to full JSON obo_ids if not.

    Args:
        adata:       AnnData object
        uberon_json: Path to UBERON JSON (used for logging only when row_ids given).
        row_ids:     Pipe-separated UBERON term IDs from the CSV row.
    """
    if not uberon_json:
        logger.info("No --uberon file provided - skipping tissue filter")
        return adata

    id_col = "tissue_ontology_term_id"
    if id_col not in adata.obs.columns:
        logger.warning(
            f"'{id_col}' not found in adata.obs — skipping tissue filter. "
            "This is expected for non-CellXGene h5ad files."
        )
        return adata

    if row_ids:
        obo_ids = set(t.strip() for t in row_ids.split("|") if t.strip())
        logger.info(f"  Tissue filter (row-specific): {sorted(obo_ids)}")
    else:
        obo_ids = load_obo_ids(uberon_json, "UBERON")

    mask     = adata.obs[id_col].isin(obo_ids)
    n_before = adata.n_obs
    adata    = adata[mask].copy()

    logger.info(f"Tissue filter: {n_before} -> {adata.n_obs} cells "
                f"({n_before - adata.n_obs} removed)")
    return adata


def filter_by_disease(adata, disease_json=None, filter_normal=False, row_ids=None):
    """Filter cells using disease_ontology_term_id.

    Uses row_ids (pipe-separated IDs from the CSV row's disease_ontology_term_id
    column) when provided. Falls back to full JSON obo_ids if not.

    Args:
        adata:         AnnData object
        disease_json:  Path to disease JSON. Skipped if None or filter_normal=False.
        filter_normal: If False, no disease filtering is applied.
        row_ids:       Pipe-separated disease ontology term IDs from the CSV row.
    """
    if not filter_normal:
        logger.info("filter_normal=False - keeping all disease states")
        return adata

    if not disease_json:
        logger.warning("filter_normal=True but no --disease file provided - skipping disease filter")
        return adata

    id_col = "disease_ontology_term_id"
    if id_col not in adata.obs.columns:
        logger.warning(
            f"'{id_col}' not found in adata.obs — skipping disease filter. "
            "This is expected for non-CellXGene h5ad files."
        )
        return adata

    if row_ids:
        obo_ids = set(t.strip() for t in row_ids.split("|") if t.strip())
        logger.info(f"  Disease filter (row-specific): {sorted(obo_ids)}")
    else:
        obo_ids  = load_obo_ids(disease_json, "disease")
    mask     = adata.obs[id_col].isin(obo_ids)
    n_before = adata.n_obs
    adata    = adata[mask].copy()

    logger.info(f"Disease filter (ontology IDs): {n_before} -> {adata.n_obs} cells "
                f"({n_before - adata.n_obs} removed)")
    return adata


def filter_by_age(adata, hsapdv_json=None, filter_normal=False, row_ids=None):
    """Filter cells using development_stage_ontology_term_id.

    Uses row_ids (pipe-separated IDs from the CSV row's
    development_stage_ontology_term_id column) when provided.
    Falls back to full JSON obo_ids if not.

    Args:
        adata:         AnnData object
        hsapdv_json:   Path to HsapDv JSON. Skipped if None or filter_normal=False.
        filter_normal: If False, no age filtering is applied.
        row_ids:       Pipe-separated HsapDv term IDs from the CSV row.
    """
    if not filter_normal:
        logger.info("filter_normal=False - skipping age filter")
        return adata

    if not hsapdv_json:
        logger.warning("filter_normal=True but no --hsapdv file provided - skipping age filter")
        return adata

    id_col = "development_stage_ontology_term_id"
    if id_col not in adata.obs.columns:
        logger.warning(
            f"'{id_col}' not found in adata.obs — skipping age filter. "
            "This is expected for non-CellXGene h5ad files."
        )
        return adata
    
    if row_ids:
        obo_ids = set(t.strip() for t in row_ids.split("|") if t.strip())
        logger.info(f"  Age filter (row-specific): {sorted(obo_ids)}")
    else:
        obo_ids    = load_obo_ids(hsapdv_json, "HsapDv")
    obo_labels = load_obo_labels(hsapdv_json)
    mask       = adata.obs[id_col].isin(obo_ids)
    n_before   = adata.n_obs

    # Log every term ID present in the data with its label, kept/excluded status and cell count
    all_counts     = adata.obs[id_col].value_counts()
    kept_terms     = [(tid, cnt) for tid, cnt in all_counts.items() if tid in obo_ids]
    excluded_terms = [(tid, cnt) for tid, cnt in all_counts.items() if tid not in obo_ids]

    logger.info(f"  Kept     ({len(kept_terms)} terms, "
                f"{sum(c for _, c in kept_terms):,} cells):")
    for term_id, count in kept_terms:
        label = obo_labels.get(term_id, "unknown label")
        logger.info(f"    KEPT     {term_id}  {label}: {count:,} cells")

    if excluded_terms:
        logger.info(f"  Excluded ({len(excluded_terms)} terms, "
                    f"{sum(c for _, c in excluded_terms):,} cells):")
        for term_id, count in excluded_terms:
            label = obo_labels.get(term_id, "unknown label — not in HsapDv JSON")
            logger.info(f"    EXCLUDED {term_id}  {label}: {count:,} cells")
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

def run_filter_adata(h5ad_path, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id,
                     filter_normal=False,
                     uberon_json=None,
                     disease_json=None,
                     hsapdv_json=None,
                     min_cluster_size=5,
                     row_uberon_ids=None,
                     row_disease_ids=None,
                     row_hsapdv_ids=None,
                     filter_obs_column=None,
                     filter_obs_value=None
):
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
        embedding:        embedding (used for output directory)
        dataset_version_id: dataset_version_id (used for output directory)
        filter_normal:    If True, apply tissue + disease + age filters
        uberon_json:      Path to UBERON JSON from resolve-uberon
        disease_json:     Path to disease JSON from resolve-disease
        hsapdv_json:      Path to HsapDv JSON from resolve-hsapdv --min-age N
        min_cluster_size: Minimum cells per cluster (default 5)
    """
    log_section("NSForest: Filter AnnData")

    prefix                = get_output_prefix( organ, first_author, journal, year, cluster_header, embedding, dataset_version_id)
    filtered_h5ad_name    = f"{prefix}_adata_filtered.h5ad"

    logger.info(f"Loading: {h5ad_path}")
    adata = load_h5ad(h5ad_path, cluster_header)
    logger.info(f"Original data: {adata.n_obs} cells, {adata.n_vars} genes, "
                f"{adata.obs[cluster_header].nunique()} clusters")

    create_stats_before_filter(adata, cluster_header, prefix, embedding)

    logger.info("\n=== Applying Filters ===")

    logger.info("\n[1/5] Tissue filter")
    adata = filter_by_tissue(adata, uberon_json, row_ids=row_uberon_ids)

    logger.info("\n[2/5] Disease filter")
    adata = filter_by_disease(adata, disease_json, filter_normal, row_ids=row_disease_ids)

    logger.info("\n[3/5] Age filter")
    adata = filter_by_age(adata, hsapdv_json, filter_normal, row_ids=row_hsapdv_ids)

    logger.info("\n[4/5] Obs column filter")
    if filter_obs_column and filter_obs_value:
        col = filter_obs_column
        val = filter_obs_value
        if col not in adata.obs.columns:
            raise ValueError(
                f"filter_obs_column '{col}' not found in adata.obs. "
                f"Available: {list(adata.obs.columns)}"
            )
        n_before = adata.n_obs
        adata = adata[adata.obs[col] == val].copy()
        logger.info(f"Obs filter adata.obs['{col}'] == '{val}': {n_before} -> {adata.n_obs} cells "
                    f"({n_before - adata.n_obs} removed)")
        if adata.n_obs == 0:
            raise ValueError(f"No cells remaining after filtering '{col}' == '{val}'")
    else:
        logger.info("No obs column filter requested — skipping")

    logger.info("\n[5/5] Min cluster size filter")
    adata = filter_by_min_cluster_size(adata, cluster_header, min_cluster_size)

    logger.info(f"\nFinal filtered data: {adata.n_obs} cells, {adata.n_vars} genes, "
                f"{adata.obs[cluster_header].nunique()} clusters")

    logger.info("\n=== Creating AFTER FILTER statistics ===")

    n_clusters = adata.obs[cluster_header].nunique()

    try:
        # Ensure X_pca exists — ns.pp.dendrogram hardcodes use_rep="X_pca"
        if "X_pca" not in adata.obsm:
            logger.info("X_pca not found in obsm — computing PCA for dendrogram")
            sc.pp.pca(adata)
        # Populate adata.uns['dendrogram_' + cluster_header] for downstream consumers.
        # SVG + CSVs are emitted by dendrogram.py (see modules/nsforest/dendrogram.nf).
        ns.pp.dendrogram(
            adata, cluster_header,
            tl_kwargs={"optimal_ordering": True}, save=False,
        )
    except ValueError as e:
        if "negative distances" in str(e):
            logger.warning(f"Optimal ordering failed — retrying without: {e}")
            ns.pp.dendrogram(
                adata, cluster_header,
                tl_kwargs={"optimal_ordering": False}, save=False,
            )
        else:
            raise
        
    # Add gene symbols into adata.var so every downstream step has them for free

    ensg_to_symbol = load_gene_mapping()
    if ensg_to_symbol:
        adata = add_gene_symbols_to_adata(adata, ensg_to_symbol)
    else:
        logger.warning("No gene mapping available — adata.var['gene_symbol'] not populated")

    adata.write_h5ad(filtered_h5ad_name)
    logger.info(f"Saved filtered h5ad: {filtered_h5ad_name}")

    logger.info("\nFilter complete!")

