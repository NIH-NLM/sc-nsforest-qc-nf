"""
Common utilities for NSForest CLI commands.
"""

import os
from pathlib import Path
import scanpy as sc
import pandas as pd
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s'
)
logger = logging.getLogger(__name__)


def create_output_dir(organ, first_author, year):
    """
    Create standardized output directory.
    
    Pattern: outputs_{organ}_{first_author}_{year}/
    """
    output_dir = f"outputs_{organ}_{first_author}_{year}"
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    return output_dir


def get_output_prefix(output_dir, cluster_header):
    """
    Get standardized output file prefix.
    
    Pattern: outputs_{organ}_{author}_{year}/{cluster_header}_
    """
    return f"{output_dir}/{cluster_header}"


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


def save_dataframe(df, filepath, formats=['csv', 'json']):
    """Save DataFrame in multiple formats."""
    for fmt in formats:
        if fmt == 'csv':
            output = f"{filepath}.csv"
            df.to_csv(output, index=False)
            logger.info(f"Saved: {output}")
        elif fmt == 'json':
            output = f"{filepath}.json"
            df.to_json(output, orient='records', indent=2)
            logger.info(f"Saved: {output}")


def log_section(title):
    """Print formatted section header."""
    logger.info("=" * 80)
    logger.info(title)
    logger.info("=" * 80)
