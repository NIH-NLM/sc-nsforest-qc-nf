"""
Common utilities for NSForest CLI commands.
"""
import csv
import scanpy as sc
import pandas as pd
import numpy as np
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

# in common_utils.py

def get_output_prefix(organ, first_author, journal, year, cluster_header, embedding="", dataset_version_id=""):
    """Build standardized output filename prefix with embedding and vid suffix for uniqueness."""
    journal_safe = journal.replace(" ", "_") if journal else "unknown"
    cluster_header_safe = cluster_header.replace(" ", "_")
    embedding_safe = embedding.replace(" ", "_") if embedding else "unknown"
    vid_suffix = f"{dataset_version_id[-6:]}" if dataset_version_id and len(dataset_version_id) >= 6 else ""
    return f"{organ}_{first_author}_{journal_safe}_{year}_{cluster_header_safe}_{embedding_safe}_{vid_suffix}"

    
def load_h5ad(h5ad_path, cluster_header):
    """Load h5ad file with validation."""
    logger.info(f"Loading h5ad: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)

    # Convert dense X to CSR sparse — large CXG datasets sometimes ship dense X
    # which inflates memory ~10x for typical sparsity. Conversion is a one-time cost.
    import scipy.sparse as sp
    if not sp.issparse(adata.X):
        if hasattr(adata.X, 'mean'):
            nnz_frac = (adata.X != 0).mean()
        else:
            nnz_frac = None

        if nnz_frac is None or nnz_frac < 0.5:
            if nnz_frac is not None:
                msg = f"Converting dense X to CSR sparse (density: {nnz_frac:.1%})"
            else:
                msg = "Converting dense X to CSR sparse"
            logger.info(msg)
            adata.X = sp.csr_matrix(adata.X)

    # Large concatenated matrices can exceed 2**31 non-zeros, which forces scipy to use int64
    # for indptr but can leave indices as int32 — a mismatch that crashes column-subset ops
    # (scipy csr_column_index1: "Output dtype not compatible with inputs"). Normalize to int64.
    if sp.issparse(adata.X) and adata.X.nnz > 2**31 - 1:
        X = adata.X.tocsr()
        if X.indices.dtype != np.int64 or X.indptr.dtype != np.int64:
            logger.info(f"Large sparse matrix (nnz={X.nnz:,} > 2^31) — promoting CSR indices to int64")
            X.indices = X.indices.astype(np.int64)
            X.indptr  = X.indptr.astype(np.int64)
            adata.X = X

    logger.info(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    if cluster_header not in adata.obs.columns:
        raise ValueError(
            f"Cluster column '{cluster_header}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    n_clusters = adata.obs[cluster_header].nunique()
    logger.info(f"Cluster header: {cluster_header} ({n_clusters} clusters)")

    return adata

def log_section(title):
    """Print formatted section header."""

    logger.info("=" * 80)
    logger.info(title)
    logger.info("=" * 80)

def save_dataframe(df, filepath, formats=['csv']):
    """Save DataFrame in specified formats."""

    for fmt in formats:
        if fmt == 'csv':
            output = f"{filepath}.csv"
            df.to_csv(output, index=False)
            logger.info(f"Saved: {output}")
        elif fmt == 'json':
            output = f"{filepath}.json"
            df.to_json(output, orient='records', indent=2)
            logger.info(f"Saved: {output}")

def setup_file_logging(name):
    fh = logging.FileHandler(f"{name}.log")
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger().addHandler(fh)


