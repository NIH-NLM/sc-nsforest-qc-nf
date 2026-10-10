"""
List the genes with a binary score above 0 in each cluster.

Reads the gene-by-cluster binary score table that prep wrote (binary_scores_ensg_{prefix}.csv)
and writes one row per cluster and gene: binary_positive_genes_{prefix}.csv with the columns
clusterName, gene_ensg, gene_symbol, binary_score, sorted by cluster and then by score
(highest first). A gene with no symbol in the filtered h5ad keeps its ENSG id as the symbol.
"""
import pandas as pd

from .common_utils import get_output_prefix, load_h5ad, log_section, logger


def run_binary_positive_genes(binary_scores_csv, filtered_h5ad, cluster_header, organ, first_author,
                              journal, year, embedding, dataset_version_id):
    log_section("NSForest: Binary Positive Genes")
    prefix = get_output_prefix(organ, first_author, journal, year, cluster_header, embedding, dataset_version_id)

    scores = pd.read_csv(binary_scores_csv, index_col=0)
    adata = load_h5ad(filtered_h5ad, cluster_header)
    if 'gene_symbol' in adata.var.columns:
        sym_map = dict(zip(adata.var_names, adata.var['gene_symbol']))
    else:
        logger.warning("adata.var['gene_symbol'] missing - the symbol column will hold the ENSG ids")
        sym_map = {}

    long = scores.stack().rename('binary_score').reset_index()
    long.columns = ['gene_ensg', 'clusterName', 'binary_score']
    long = long[long['binary_score'] > 0].copy()
    long['gene_symbol'] = [sym_map.get(g, g) for g in long['gene_ensg']]
    long = long[['clusterName', 'gene_ensg', 'gene_symbol', 'binary_score']]
    long = long.sort_values(['clusterName', 'binary_score'], ascending=[True, False])

    out = f"binary_positive_genes_{prefix}.csv"
    long.to_csv(out, index=False)
    logger.info(f"Saved: {out} ({len(long)} rows, {long['clusterName'].nunique()} clusters)")
