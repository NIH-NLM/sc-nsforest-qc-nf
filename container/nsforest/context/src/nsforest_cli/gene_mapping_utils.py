"""
Gene mapping utilities for converting Ensembl IDs to gene symbols.

Reads gene mapping from NIH-NLM cell-kn GitHub repository.
"""

import pandas as pd
import requests
from io import StringIO

from .common_utils import logger

GENE_MAPPING_URL = "https://raw.githubusercontent.com/NIH-NLM/cell-kn/main/data/biomart/gene_mapping.csv"


def load_gene_mapping(url=None):
    """
    Load gene mapping CSV from GitHub.
    
    Parameters
    ----------
    url : str, optional
        URL to gene mapping CSV. Defaults to cell-kn repository.
        
    Returns
    -------
    dict
        Dictionary mapping ensembl_gene_id → external_gene_name
    """
    if url is None:
        url = GENE_MAPPING_URL
    
    logger.info(f"Loading gene mapping from: {url}")
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        df = pd.read_csv(StringIO(response.text))
        
        # Create mapping dict
        ensg_to_symbol = dict(zip(df['ensembl_gene_id'], df['external_gene_name']))
        
        logger.info(f"Loaded {len(ensg_to_symbol)} gene mappings")
        
        return ensg_to_symbol
        
    except Exception as e:
        logger.error(f"Failed to load gene mapping: {e}")
        logger.warning("Continuing without gene mapping - will use Ensembl IDs")
        return {}


def map_markers_to_symbols(results_df, ensg_to_symbol):
    """
    Map NSForest marker Ensembl IDs to gene symbols.
    
    Parameters
    ----------
    results_df : pandas.DataFrame
        NSForest results with 'NSForest_markers' column
    ensg_to_symbol : dict
        Mapping of ensembl_gene_id → external_gene_name
        
    Returns
    -------
    pandas.DataFrame
        Results with added 'gene_names' column
    dict
        markers_dict: cluster name → gene symbols list
    """
    logger.info("Mapping NSForest markers to gene symbols...")
    
    results_df['gene_names'] = [
        [ensg_to_symbol.get(gene, gene) for gene in markers]
        for markers in results_df['NSForest_markers']
    ]
    
    # Create markers_dict for plotting
    markers_dict = dict(zip(results_df["clusterName"], results_df["gene_names"]))
    
    logger.info(f"Mapped markers for {len(results_df)} clusters")
    
    return results_df, markers_dict


def add_gene_symbols_to_adata(adata, ensg_to_symbol):
    """
    Add gene_symbol column to adata.var for plotting.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    ensg_to_symbol : dict
        Mapping of ensembl_gene_id → external_gene_name
        
    Returns
    -------
    AnnData
        adata with 'gene_symbol' column in .var
    """
    logger.info("Adding gene symbols to adata.var...")
    
    adata.var['gene_symbol'] = [
        ensg_to_symbol.get(gene, gene) for gene in adata.var_names
    ]
    
    logger.info("Gene symbols added to adata.var['gene_symbol']")
    
    return adata
