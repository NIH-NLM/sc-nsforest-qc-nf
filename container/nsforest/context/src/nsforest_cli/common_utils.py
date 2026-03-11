"""
Common utilities for NSForest CLI commands.
"""

import scanpy as sc
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

# in common_utils.py

def get_output_prefix(organ, first_author, year, cluster_header):
    output_prefix = f"{organ}_{first_author}_{year}_{cluster_header}"
    return output_prefix 

    
def load_h5ad(h5ad_path, cluster_header):
    """Load h5ad file with validation."""

    logger.info(f"Loading h5ad: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
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


