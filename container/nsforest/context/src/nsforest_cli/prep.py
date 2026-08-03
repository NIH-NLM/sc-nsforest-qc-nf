"""
Compute median expression AND binary scores per cluster in a single pass.

Loads the h5ad once, runs ns.pp.prep_medians() then ns.pp.prep_binary_scores()
on the same in-memory object (matching DEMO_NS-Forest_workflow), and writes both
sets of outputs. Replaces the separate prep-medians / prep-binary-scores pipeline
steps so the (expensive, densifying) median computation happens exactly once.

Saves:
  medians_ensg_{prefix}.csv/.pkl        (+ medians_symbols_{prefix}.csv/.pkl)
  binary_scores_ensg_{prefix}.csv/.pkl  (+ binary_scores_symbols_{prefix}.csv/.pkl)
"""
import nsforest as ns
import pandas as pd

from .common_utils import get_output_prefix, load_h5ad, log_section, logger


def _varm_to_df(adata_prep, key):
    """Rebuild a gene-by-cluster DataFrame with an explicit ENSG index.

    Some anndata versions return a varm DataFrame whose index assignment doesn't
    persist; build a fresh one from var_names (matches run_nsforest's expectations).
    """
    raw = adata_prep.varm[key]
    values = raw.values if hasattr(raw, 'values') else raw
    columns = raw.columns if hasattr(raw, 'columns') else None
    return pd.DataFrame(
        values,
        index=pd.Index(list(adata_prep.var_names), name='gene'),
        columns=columns,
    )


def _write_pair(df, sym_map, prefix, kind):
    """Write ensg + symbols csv/pkl for one flavor ('medians' or 'binary_scores')."""
    df.to_csv(f"{kind}_ensg_{prefix}.csv")
    df.to_pickle(f"{kind}_ensg_{prefix}.pkl")
    logger.info(f"Saved: {kind}_ensg_{prefix}.csv / .pkl")

    if sym_map is not None:
        df_sym = df.rename(index=lambda g: sym_map.get(g, g))
        df_sym.to_csv(f"{kind}_symbols_{prefix}.csv")
        df_sym.to_pickle(f"{kind}_symbols_{prefix}.pkl")
        logger.info(f"Saved: {kind}_symbols_{prefix}.csv / .pkl")
    else:
        logger.warning(f"gene_symbol unavailable — skipping {kind}_symbols outputs")


def run_prep(h5ad_path, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id):
    """Compute medians and binary scores in one pass and write both output sets."""
    log_section("NSForest: Prep (medians + binary scores, single pass)")
    prefix = get_output_prefix(organ, first_author, journal, year, cluster_header, embedding, dataset_version_id)

    adata = load_h5ad(h5ad_path, cluster_header)

    # Capture the ENSG->symbol map up front, from the freshly loaded adata, so the
    # symbol outputs never depend on `adata` surviving the prep calls below.
    if 'gene_symbol' in adata.var.columns:
        sym_map = dict(zip(adata.var_names, adata.var['gene_symbol']))
    else:
        sym_map = None
        logger.warning("adata.var['gene_symbol'] missing — symbol outputs will be skipped")

    # One densifying pass. No .copy(): prep_medians subsets positive genes into a NEW
    # object and only adds to varm — it never mutates adata.X — so a copy just doubles peak RAM.
    logger.info("Running ns.pp.prep_medians()...")
    adata_prep = ns.pp.prep_medians(adata, cluster_header)
    df_medians = _varm_to_df(adata_prep, 'medians_' + cluster_header)
    logger.info(f"Medians shape: {df_medians.shape}, index dtype: {df_medians.index.dtype}")

    # Binary scores read ONLY varm['medians_...'] (never adata.X) — no cell densification here.
    logger.info("Running ns.pp.prep_binary_scores()...")
    adata_prep = ns.pp.prep_binary_scores(adata_prep, cluster_header)
    df_binary = _varm_to_df(adata_prep, 'binary_scores_' + cluster_header)
    logger.info(f"Binary scores shape: {df_binary.shape}")

    _write_pair(df_medians, sym_map, prefix, 'medians')
    _write_pair(df_binary, sym_map, prefix, 'binary_scores')

    logger.info("Prep complete!")
