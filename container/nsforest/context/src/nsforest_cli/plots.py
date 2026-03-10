"""
Create NSForest visualization plots.

Corresponds to DEMO_NS-Forest_workflow.py: Section 4

Saves boxplots (html), scatter plots (svg), and expression plots (svg).
"""

import matplotlib
matplotlib.use("Agg")

import ast
import pandas as pd
import nsforest as ns

from .common_utils import (
    load_h5ad,
    log_section,
    logger
)
from .gene_mapping_utils import (
    load_gene_mapping,
    map_markers_to_symbols,
    add_gene_symbols_to_adata
)


def run_plots(h5ad_path, results_csv, cluster_header, organ, first_author, year):
    """
    Create NSForest visualization plots with gene symbol mapping.
    """
    log_section("NSForest: Plotting")

    cluster_header_safe = cluster_header.replace(" ", "_")
    prefix = f"{organ}_{first_author}_{year}_{cluster_header_safe}"

    # Load results
    logger.info(f"Loading results: {results_csv}")
    results = pd.read_csv(results_csv)
    results['NSForest_markers'] = results['NSForest_markers'].apply(ast.literal_eval)

    # Gene symbol mapping
    ensg_to_symbol = load_gene_mapping()
    results, markers_dict = map_markers_to_symbols(results, ensg_to_symbol)

    # Boxplots (html)
    for metric in ['f_score', 'precision', 'recall', 'onTarget']:
        ns.pl.boxplot(results, metric, save="html", output_folder="", outputfilename_prefix=prefix)
    logger.info("Boxplots saved.")

    # Scatter plots (svg)
    for metric in ['f_score', 'precision', 'recall', 'onTarget']:
        ns.pl.scatter_w_clusterSize(results, metric, save=True, output_folder="", outputfilename_prefix=prefix)
    logger.info("Scatter plots saved.")

    # Load adata for expression plots
    logger.info(f"Loading h5ad: {h5ad_path}")
    adata = load_h5ad(h5ad_path, cluster_header)
    adata = add_gene_symbols_to_adata(adata, ensg_to_symbol)

    # Dotplot
    ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                  gene_symbols='gene_symbol', save="svg", output_folder="",
                  outputfilename_suffix=prefix)
    ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                  gene_symbols='gene_symbol', standard_scale='var', save="svg",
                  output_folder="", outputfilename_suffix=prefix + "_scaled")

    # Stacked violin
    ns.pl.stackedviolin(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                        gene_symbols='gene_symbol', save="svg", output_folder="",
                        outputfilename_suffix=prefix)
    ns.pl.stackedviolin(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                        gene_symbols='gene_symbol', standard_scale='var', save="svg",
                        output_folder="", outputfilename_suffix=prefix + "_scaled")

    # Matrix plot
    ns.pl.matrixplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                     gene_symbols='gene_symbol', save="svg", output_folder="",
                     outputfilename_suffix=prefix)
    ns.pl.matrixplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                     gene_symbols='gene_symbol', standard_scale='var', save="svg",
                     output_folder="", outputfilename_suffix=prefix + "_scaled")

    logger.info("Plotting complete!")
