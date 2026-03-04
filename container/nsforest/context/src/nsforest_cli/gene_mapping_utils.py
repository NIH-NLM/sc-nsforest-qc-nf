"""
Gene mapping utilities for converting Ensembl IDs to gene symbols.

Reads gene mapping from baked-in local file, falls back to GitHub if not found.
"""

import os
import pandas as pd
import requests
from io import StringIO

from .common_utils import logger

LOCAL_GENE_MAPPING  = "/app/gene_mapping.csv"
GENE_MAPPING_URL    = "https://raw.githubusercontent.com/NIH-NLM/cell-kn/main/data/prod/biomart/gene_mapping.csv"


def load_gene_mapping(url=None):
    """
    Load gene mapping CSV — local file first, GitHub fallback.

    Returns
    -------
    dict
        Dictionary mapping ensembl_gene_id → external_gene_name
    """
    if os.path.exists(LOCAL_GENE_MAPPING):
        logger.info(f"Loading gene mapping from local file: {LOCAL_GENE_MAPPING}")
        try:
            df = pd.read_csv(LOCAL_GENE_MAPPING)
            ensg_to_symbol = dict(zip(df['ensembl_gene_id'], df['external_gene_name']))
            logger.info(f"Loaded {len(ensg_to_symbol)} gene mappings")
            return ensg_to_symbol
        except Exception as e:
            logger.warning(f"Failed to read local gene mapping: {e} — trying network")

    # Fallback to network
    if url is None:
        url = GENE_MAPPING_URL

    logger.info(f"Loading gene mapping from network: {url}")
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        df = pd.read_csv(StringIO(response.text))
        ensg_to_symbol = dict(zip(df['ensembl_gene_id'], df['external_gene_name']))
        logger.info(f"Loaded {len(ensg_to_symbol)} gene mappings")
        return ensg_to_symbol
    except Exception as e:
        logger.error(f"Failed to load gene mapping from network: {e}")
        logger.warning("Continuing without gene mapping — will use Ensembl IDs")
        return {}


def map_markers_to_symbols(results_df, ensg_to_symbol):
    """Map NSForest marker Ensembl IDs to gene symbols."""
    logger.info("Mapping NSForest markers to gene symbols...")

    results_df['gene_names'] = [
        [ensg_to_symbol.get(gene, gene) for gene in markers]
        for markers in results_df['NSForest_markers']
    ]

    markers_dict = dict(zip(results_df["clusterName"], results_df["gene_names"]))
    logger.info(f"Mapped markers for {len(results_df)} clusters")
    return results_df, markers_dict


def add_gene_symbols_to_adata(adata, ensg_to_symbol):
    """Add gene_symbol column to adata.var for plotting."""
    logger.info("Adding gene symbols to adata.var...")
    adata.var['gene_symbol'] = [
        ensg_to_symbol.get(gene, gene) for gene in adata.var_names
    ]
    logger.info("Gene symbols added to adata.var['gene_symbol']")
    return adata
