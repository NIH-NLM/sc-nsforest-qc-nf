"""
Create NSForest visualization plots.

Generates boxplots, scatter plots, and expression plots with gene symbol mapping.
"""

import ast
import pandas as pd
import nsforest as ns

from .common_utils import (
    create_output_dir,
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
    
    Reads:
    - h5ad file (original adata for expression plots)
    - results.csv (NSForest results)
    
    Creates:
    Boxplots (HTML):
    - boxplot_f_score_{cluster_header}.html
    - boxplot_precision_{cluster_header}.html
    - boxplot_recall_{cluster_header}.html
    - boxplot_onTarget_{cluster_header}.html
    
    Scatter plots with cluster size (SVG):
    - scatter_w_clusterSize_f_score_{cluster_header}.svg
    - scatter_w_clusterSize_precision_{cluster_header}.svg
    - scatter_w_clusterSize_recall_{cluster_header}.svg
    - scatter_w_clusterSize_onTarget_{cluster_header}.svg
    
    Expression plots (SVG):
    - dotplot_{cluster_header}.svg
    - dotplot_{cluster_header}_scaled.svg
    - stackedviolin_{cluster_header}.svg
    - stackedviolin_{cluster_header}_scaled.svg
    - matrixplot_{cluster_header}.svg
    - matrixplot_{cluster_header}_scaled.svg
    """
    log_section("NSForest: Plotting")
    
    output_folder = create_output_dir(organ, first_author, year) + "/"
    outputfilename_prefix = cluster_header
    
    # Load results
    logger.info(f"Loading results: {results_csv}")
    results = pd.read_csv(results_csv)
    logger.info(f"Results shape: {results.shape}")
    
    # CRITICAL: Parse NSForest_markers from string to list
    logger.info("Parsing NSForest_markers from CSV string format...")
    results['NSForest_markers'] = results['NSForest_markers'].apply(ast.literal_eval)

    # Load gene mapping
    ensg_to_symbol = load_gene_mapping()
    
    # Map markers to gene symbols
    results, markers_dict = map_markers_to_symbols(results, ensg_to_symbol)
    
    # Boxplots (HTML)
    metrics = ['f_score', 'precision', 'recall', 'onTarget']
    
    logger.info("Creating boxplots...")
    for metric in metrics:
        logger.info(f"  Boxplot: {metric}")
        ns.pl.boxplot(
            results, 
            metric, 
            save="html", 
            output_folder=output_folder, 
            outputfilename_prefix=outputfilename_prefix
        )
    
    # Scatter plots with cluster size (SVG)
    logger.info("Creating scatter plots with cluster size...")
    for metric in metrics:
        logger.info(f"  Scatter: {metric}")
        ns.pl.scatter_w_clusterSize(
            results, 
            metric, 
            save=True, 
            output_folder=output_folder, 
            outputfilename_prefix=outputfilename_prefix
        )
    
    # Expression plots
    logger.info(f"Loading h5ad for expression plots: {h5ad_path}")
    adata = load_h5ad(h5ad_path, cluster_header)
    
    # Add gene symbols to adata
    adata = add_gene_symbols_to_adata(adata, ensg_to_symbol)
    
    # Dotplot
    logger.info("Creating dotplot...")
    ns.pl.dotplot(
        adata, 
        markers_dict, 
        cluster_header, 
        dendrogram=True, 
        use_raw=False,
        gene_symbols='gene_symbol',
        save="svg", 
        output_folder=output_folder, 
        outputfilename_suffix=outputfilename_prefix
    )
    
    logger.info("Creating dotplot (scaled)...")
    ns.pl.dotplot(
        adata, 
        markers_dict, 
        cluster_header, 
        dendrogram=True, 
        use_raw=False,
        gene_symbols='gene_symbol',
        standard_scale='var',
        save="svg", 
        output_folder=output_folder, 
        outputfilename_suffix=outputfilename_prefix + "_scaled"
    )
    
    # Stacked violin
    logger.info("Creating stacked violin plot...")
    ns.pl.stackedviolin(
        adata, 
        markers_dict, 
        cluster_header, 
        dendrogram=True, 
        use_raw=False,
        gene_symbols='gene_symbol',
        save="svg", 
        output_folder=output_folder, 
        outputfilename_suffix=outputfilename_prefix
    )
    
    logger.info("Creating stacked violin plot (scaled)...")
    ns.pl.stackedviolin(
        adata, 
        markers_dict, 
        cluster_header, 
        dendrogram=True, 
        use_raw=False,
        gene_symbols='gene_symbol',
        standard_scale='var',
        save="svg", 
        output_folder=output_folder, 
        outputfilename_suffix=outputfilename_prefix + "_scaled"
    )
    
    # Matrix plot
    logger.info("Creating matrix plot...")
    ns.pl.matrixplot(
        adata, 
        markers_dict, 
        cluster_header, 
        dendrogram=True, 
        use_raw=False,
        gene_symbols='gene_symbol',
        save="svg", 
        output_folder=output_folder, 
        outputfilename_suffix=outputfilename_prefix
    )
    
    logger.info("Creating matrix plot (scaled)...")
    ns.pl.matrixplot(
        adata, 
        markers_dict, 
        cluster_header, 
        dendrogram=True, 
        use_raw=False,
        gene_symbols='gene_symbol',
        standard_scale='var',
        save="svg", 
        output_folder=output_folder, 
        outputfilename_suffix=outputfilename_prefix + "_scaled"
    )
    
    logger.info("Plotting complete!")
