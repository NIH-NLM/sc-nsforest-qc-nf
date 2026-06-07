"""
Compute median expression per cluster.

Corresponds to DEMO_NS-Forest_workflow.py: Section 3 prep
Uses ns.pp.prep_medians() to filter positive genes and compute medians.

Saves:
  medians_{organ}_{first_author}_{journal}_{year}_{cluster_header}_{embedding}_{vid}.csv
  medians_symbols_{organ}_{first_author}_{journal}_{year}_{cluster_header}_{embedding}_{vid}.csv
  medians_{organ}_{first_author}_{journal}_{year}_{cluster_header}_{embedding}_{vid}.pkl
  medians_symbols_{organ}_{first_author}_{journal}_{year}_{cluster_header}_{embedding}_{vid}.pkl
"""
import csv
import nsforest as ns
import pandas as pd

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
    raw = adata_prep.varm['medians_' + cluster_header]

    # Bypass anndata varm wrapping — build a fresh DataFrame with explicit ENSG index.
    # Some anndata versions return a DataFrame whose index assignment doesn't persist.
    values = raw.values if hasattr(raw, 'values') else raw
    columns = raw.columns if hasattr(raw, 'columns') else None
    df_medians = pd.DataFrame(
        values,
        index=pd.Index(list(adata_prep.var_names), name='gene'),
        columns=columns,
    )

    logger.info(f"Medians shape: {df_medians.shape}, index dtype: {df_medians.index.dtype}")
    logger.info(f"First 3 gene IDs: {list(df_medians.index[:3])}")
    
    if 'gene_symbol' in adata.var.columns:
        sym_map = dict(zip(adata.var_names, adata.var['gene_symbol']))
        df_medians_symbols = df_medians.rename(index=lambda g: sym_map.get(g, g))
        df_medians_symbols.to_csv(f"medians_symbols_{prefix}.csv")
        df_medians_symbols.to_pickle(f"medians_symbols_{prefix}.pkl")
        logger.info(f"Saved: medians_symbols_{prefix}.csv")
        logger.info(f"Saved: medians_symbols_{prefix}.pkl")
    else:
        logger.warning("adata.var['gene_symbol'] missing — skipping medians _symbols outputs")
    
    df_medians.to_csv(f"medians_ensg_{prefix}.csv")
    df_medians.to_pickle(f"medians_ensg_{prefix}.pkl")
    logger.info(f"Saved: medians_ensg_{prefix}.csv")
    logger.info(f"Saved: medians_ensg_{prefix}.pkl")

    logger.info("Prep medians complete!")
