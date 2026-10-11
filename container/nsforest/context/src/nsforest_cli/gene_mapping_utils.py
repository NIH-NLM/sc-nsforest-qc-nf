"""
Gene names for Ensembl ids.

The names come from the table of one GENCODE release, gene_mapping_gencode_v<N>.csv, which is
in the image (see cellxgene-harvester, build-gene-mapping). Every gene has a name: GENCODE gives
the gene id as the name of a gene that has no symbol, and a gene that is not in the table keeps
its own id. So a column of gene names never has a blank and never drops a gene.
"""

import glob
import json
import os
import re

import pandas as pd

from .common_utils import logger

GENE_MAPPING_FOLDER = "/app"
GENE_MAPPING_PATTERN = "gene_mapping_gencode_v*.csv"


def mapping_file(folder=GENE_MAPPING_FOLDER):
    """The table of the image: gene_mapping_gencode_v<N>.csv, the highest release if there are
    several. Raises FileNotFoundError if there is none: the symbol files must not be written
    without it."""
    found = glob.glob(os.path.join(folder, GENE_MAPPING_PATTERN))
    if not found:
        raise FileNotFoundError(f"no {GENE_MAPPING_PATTERN} in {folder}: the image has no gene name table")
    return max(found, key=release_of)


def release_of(path):
    """The GENCODE release in the file name: gene_mapping_gencode_v48.csv gives 48."""
    return int(re.search(r"_v(\d+)\.csv$", path).group(1))


def load_gene_mapping(folder=GENE_MAPPING_FOLDER):
    """The dict Ensembl id to gene name of the GENCODE release table in the image."""
    path = mapping_file(folder)
    df = pd.read_csv(path, usecols=["ensembl_gene_id", "gene_name"])
    ensg_to_symbol = dict(zip(df["ensembl_gene_id"], df["gene_name"]))
    logger.info(f"Loaded {len(ensg_to_symbol)} gene names from {path} (GENCODE v{release_of(path)})")
    return ensg_to_symbol


def describe_mapping(folder=GENE_MAPPING_FOLDER):
    """What the ETL JSON records about the table in use: file, GENCODE release, rows, SHA-256 and
    the provenance file that belongs to it (source, checksums, date)."""
    import hashlib
    path = mapping_file(folder)
    provenance_path = path[:-len(".csv")] + ".provenance.json"
    provenance = json.load(open(provenance_path, encoding="utf-8")) if os.path.exists(provenance_path) else None
    with open(path, "rb") as f:
        sha256 = hashlib.sha256(f.read()).hexdigest()
    return {"file": os.path.basename(path), "gencode_release": release_of(path),
            "rows": len(pd.read_csv(path, usecols=["ensembl_gene_id"])), "sha256": sha256,
            "provenance": provenance}


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
