"""
Extract per-cluster cell-ontology mapping table.

Emits a CSV with exactly four columns:
    cluster_name, skos, manual_mapped_cid, cell_ontology_id

- cluster_name       : unique value from the cluster-header column in obs.
- skos               : restricted vocabulary -- only "exact" or "related"
                       (blank by default; fill in manually or pass --default-skos).
- manual_mapped_cid  : blank, to be filled in manually downstream.
- cell_ontology_id   : value from the cell-ontology obs column
                       (default: cell_type_ontology_term_id).

One row per (cluster, cell_ontology_id) pair so multi-mapped clusters are not
silently collapsed.
"""

import pandas as pd

from .common_utils import (
    get_output_prefix,
    load_h5ad,
    log_section,
    logger,
    save_dataframe,
    setup_file_logging,
)

ALLOWED_SKOS = {"", "exact", "related"}


def build_cluster_cid_mapping(adata, cluster_header, cid_column, default_skos=""):
    """Build the 4-column cluster -> cell_ontology_id mapping dataframe."""

    if default_skos not in ALLOWED_SKOS:
        raise ValueError(
            f"default_skos must be one of {sorted(ALLOWED_SKOS)}, got {default_skos!r}"
        )

    if cid_column not in adata.obs.columns:
        raise ValueError(
            f"cell-ontology column '{cid_column}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    pairs = (
        adata.obs[[cluster_header, cid_column]]
        .astype(str)
        .drop_duplicates()
        .sort_values([cluster_header, cid_column])
        .reset_index(drop=True)
    )

    logger.info(
        f"Found {pairs[cluster_header].nunique()} unique clusters, "
        f"{len(pairs)} (cluster, cid) pairs"
    )

    return pd.DataFrame(
        {
            "cluster_name": pairs[cluster_header],
            "skos": default_skos,
            "manual_mapped_cid": "",
            "cell_ontology_id": pairs[cid_column],
        }
    )


def run_cluster_cid_mapping(
    h5ad_path,
    cluster_header,
    organ,
    first_author,
    journal,
    year,
    embedding,
    dataset_version_id,
    cid_column="cell_type_ontology_term_id",
    default_skos="",
):
    """Main entry point: load h5ad, build mapping, save CSV."""

    setup_file_logging("cluster_cid_mapping")
    log_section("NSForest: Cluster -> Cell Ontology Mapping")

    output_prefix = get_output_prefix(
        organ, first_author, journal, year, cluster_header, embedding, dataset_version_id
    )

    adata = load_h5ad(h5ad_path, cluster_header)

    mapping_df = build_cluster_cid_mapping(
        adata,
        cluster_header=cluster_header,
        cid_column=cid_column,
        default_skos=default_skos,
    )

    out_path = f"{output_prefix}_cluster_cid_mapping"
    save_dataframe(mapping_df, out_path, formats=["csv"])

    logger.info("Cluster -> cell ontology mapping complete!")
