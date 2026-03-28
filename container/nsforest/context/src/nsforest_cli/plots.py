"""
Create NSForest visualization plots.

Corresponds to DEMO_NS-Forest_workflow.py: Section 4

Saves boxplots (html), scatter plots (svg), and expression plots (svg).
"""

import matplotlib
matplotlib.use("Agg")

import ast
import glob
import os
import pandas as pd
import nsforest as ns

from .common_utils import (
    get_output_prefix,
    load_h5ad,
    log_section,
    logger
)
from .gene_mapping_utils import (
    load_gene_mapping,
    map_markers_to_symbols,
    add_gene_symbols_to_adata
)


def run_plots(h5ad_path, results_csv, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id):
    """
    Create NSForest visualization plots with gene symbol mapping.
    """
    log_section("NSForest: Plotting")

    prefix = get_output_prefix( organ, first_author, journal, year, cluster_header, embedding, dataset_version_id )

    # Load results
    logger.info(f"Loading results: {results_csv}")

    results = pd.read_csv(results_csv)
    results = results.dropna(subset=['NSForest_markers'])
    results = results[results['NSForest_markers'].str.strip() != '[]']

    if results.empty:
        logger.warning("No valid marker results — skipping all plots")
        return

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

    # Rename expression plot SVGs: move prefix from suffix to prefix position
    # NSForest library creates: {plot_type}_{prefix}[_scaled].svg
    # We want:                  {prefix}_{plot_type}[_scaled].svg
    for svg in glob.glob("*.svg"):
        if not svg.startswith(prefix):
            base = svg
            base = base.replace(f"_{prefix}_scaled.svg", "_scaled.svg")
            base = base.replace(f"_{prefix}.svg", ".svg")
            new_name = f"{prefix}_{base}"
            if new_name != svg and not os.path.exists(new_name):
                os.rename(svg, new_name)
                logger.info(f"Renamed: {svg} -> {new_name}")

    logger.info("Plotting complete!")
