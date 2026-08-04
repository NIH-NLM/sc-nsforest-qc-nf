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
import scanpy as sc
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


def run_plots(h5ad_path, results_csv, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id,
              max_cells_per_cluster=0, seed=42):
    """
    Create NSForest visualization plots with gene symbol mapping.
    """
    log_section("NSForest: Plotting")
    sc.settings.figdir = "."

    prefix = get_output_prefix( organ, first_author, journal, year, cluster_header, embedding, dataset_version_id )

    # Load results
    logger.info(f"Loading results: {results_csv}")

    results = pd.read_csv(results_csv)
    results = results.dropna(subset=['NSForest_markers', 'clusterName'])
    results = results[results['NSForest_markers'].str.strip() != '[]']

    if results.empty:
        logger.warning("No valid marker results — skipping all plots")
        return

    results['NSForest_markers'] = results['NSForest_markers'].apply(ast.literal_eval)
    
    # Gene symbol mapping
    ensg_to_symbol = load_gene_mapping()
    results, markers_dict = map_markers_to_symbols(results, ensg_to_symbol)

    # Boxplots — html (interactive) + svg (static for publication)
    for metric in ['f_score', 'precision', 'recall', 'onTarget']:
        ns.pl.boxplot(results, metric, save=True, output_folder="", outputfilename_prefix=prefix)
        os.rename(f"{prefix}_boxplot_{metric}.html", f"boxplot_{prefix}_{metric}.html")
        ns.pl.boxplot(results, metric, save="svg",  output_folder="", outputfilename_prefix=prefix)
        os.rename(f"{prefix}_boxplot_{metric}.svg",  f"boxplot_{prefix}_{metric}.svg")
    logger.info("Boxplots saved.")

    # Scatter plots — html (interactive) + svg (static for publication)
    for metric in ['f_score', 'precision', 'recall', 'onTarget']:
        ns.pl.scatter_w_clusterSize(results, metric, save=True,  output_folder="", outputfilename_prefix=prefix)
        os.rename(f"{prefix}_scatter_{metric}.html", f"scatter_{prefix}_{metric}.html")
        ns.pl.scatter_w_clusterSize(results, metric, save="svg", output_folder="", outputfilename_prefix=prefix)
        os.rename(f"{prefix}_scatter_{metric}.svg",  f"scatter_{prefix}_{metric}.svg")
    logger.info("Scatter plots saved.")
    
    # Load adata for expression plots
    logger.info(f"Loading h5ad: {h5ad_path}")
    adata = load_h5ad(h5ad_path, cluster_header)
    adata = add_gene_symbols_to_adata(adata, ensg_to_symbol)

    # some cells in adata.obs[cluster_header] have NaN (float) instead of a string label.
    # The dendrogram reorder will fail because it cant join floats as strings...
    adata = adata[adata.obs[cluster_header].notna()].copy()
    adata.obs[cluster_header] = adata.obs[cluster_header].astype(str).astype("category")

    # Cap cells per cluster for the expression plots — dotplot/violin/matrix aggregate per cluster,
    # so a few thousand cells/cluster is plenty, while the full multi-million-cell matrix OOMs
    # plotting. Stratified + seeded; clusters smaller than the cap keep all their cells.
    if max_cells_per_cluster and max_cells_per_cluster > 0:
        import numpy as np
        rng = np.random.default_rng(seed)
        groups = adata.obs.groupby(cluster_header, observed=True).indices
        keep = []
        for cl, idx in groups.items():
            idx = np.asarray(idx)
            if len(idx) > max_cells_per_cluster:
                idx = rng.choice(idx, size=max_cells_per_cluster, replace=False)
            keep.append(idx)
        keep_idx = np.sort(np.concatenate(keep))
        n_before = adata.n_obs
        adata = adata[keep_idx].copy()
        logger.info(f"Subsampled cells for plotting: {n_before} -> {adata.n_obs} "
                    f"(max {max_cells_per_cluster}/cluster, seed={seed})")

    # Dotplot
    ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                  gene_symbols='gene_symbol', save="svg", output_folder="",
                  outputfilename_suffix=prefix)
    ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                  gene_symbols='gene_symbol', standard_scale='var', save="svg",
                  output_folder="", outputfilename_suffix="_scaled" + prefix)

    # Stacked violin
    ns.pl.stackedviolin(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                        gene_symbols='gene_symbol', save="svg", output_folder="",
                        outputfilename_suffix=prefix)
    ns.pl.stackedviolin(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                        gene_symbols='gene_symbol', standard_scale='var', save="svg",
                        output_folder="", outputfilename_suffix="_scaled" + prefix)

    # Matrix plot
    ns.pl.matrixplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                     gene_symbols='gene_symbol', save="svg", output_folder="",
                     outputfilename_suffix=prefix)
    ns.pl.matrixplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                     gene_symbols='gene_symbol', standard_scale='var', save="svg",
                     output_folder="", outputfilename_suffix="_scaled" + prefix)


    logger.info("Plotting complete!")
