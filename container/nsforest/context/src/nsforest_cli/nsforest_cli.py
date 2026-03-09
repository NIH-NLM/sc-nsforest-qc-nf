#!/usr/bin/env python
"""
NSForest CLI - Marker gene discovery for single-cell data

Workflow matches DEMO_NS-forest_workflow.ipynb
"""

import typer
from pathlib import Path

app = typer.Typer(add_completion=False, help="NSForest marker discovery workflow")

@app.command("dendrogram")
def dendrogram_command(
    h5ad_path: Path = typer.Option(..., help="Path to input h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """Generate hierarchical clustering dendrogram.
       Corresponds to DEMO notebook: Section 1 - Data loading and dendrogram
    """
    from dendrogram import run_dendrogram
    run_dendrogram(h5ad_path, cluster_header, organ, first_author, year)


@app.command("cluster-stats")
def cluster_stats_command(
    h5ad_path: Path = typer.Option(..., help="Path to input h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """
    Compute cluster statistics (cell counts, percentages).

    Corresponds to DEMO notebook: Section 2 - Cluster statistics
    """
    from cluster_stats import run_cluster_stats
    run_cluster_stats(h5ad_path, cluster_header, organ, first_author, year)


@app.command("prep-medians-binary-scores")
def prep_medians_binary_scores_command(
    h5ad_path: Path = typer.Option(..., help="Path to input h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """
    Compute median expression and binary scores per cluster.

    Runs ns.pp.prep_medians() then ns.pp.prep_binary_scores() on the same
    adata_prep object, matching DEMO notebook Section 3 exactly.

    Saves: adata_prep.h5ad, {cluster_header}_medians.csv, {cluster_header}_binary_scores.csv
    """
    from prep_medians_binary_scores import run_prep_medians_binary_scores
    run_prep_medians_binary_scores(h5ad_path, cluster_header, organ, first_author, year)


@app.command("plot-histograms")
def plot_histograms_command(
    binary_scores_path: Path = typer.Option(..., help="Path to binary scores CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """
    Plot histograms of binary score distributions.

    Corresponds to DEMO notebook: Section 5 - Plot histograms
    """
    from commands.plot_histograms import run_plot_histograms
    run_plot_histograms(binary_scores_path, cluster_header, organ, first_author, year)


@app.command("run-nsforest")
def run_nsforest_command(
    binary_scores_path: Path = typer.Option(..., help="Path to binary scores CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
    n_trees: int = typer.Option(1000, help="Number of trees in random forest"),
    n_jobs: int = typer.Option(4, help="Number of parallel jobs"),
    cluster_list: str = typer.Option(None, help="Comma-separated list of clusters (for parallelization)"),
):
    """
    Run NSForest random forest marker discovery.

    Corresponds to DEMO notebook: Section 6 - Run NSForest
    Can be parallelized across clusters using --cluster-list
    """
    from commands.run_nsforest import run_nsforest
    clusters = cluster_list.split(',') if cluster_list else None
    run_nsforest(binary_scores_path, cluster_header, organ, first_author, year,
                 n_trees, n_jobs, clusters)


@app.command("plot-boxplots")
def plot_boxplots_command(
    results_path: Path = typer.Option(..., help="Path to NSForest results CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """
    Plot boxplots of feature importance.

    Corresponds to DEMO notebook: Section 7 - Plot boxplots
    """
    from commands.plot_boxplots import run_plot_boxplots
    run_plot_boxplots(results_path, cluster_header, organ, first_author, year)


@app.command("plot-scatter")
def plot_scatter_command(
    results_path: Path = typer.Option(..., help="Path to NSForest results CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """
    Plot scatter of cluster separation metrics.

    Corresponds to DEMO notebook: Section 8 - Plot scatter
    """
    from commands.plot_scatter import run_plot_scatter
    run_plot_scatter(results_path, cluster_header, organ, first_author, year)


@app.command("gene-mapping")
def gene_mapping_command(
    h5ad_path: Path = typer.Option(..., help="Path to input h5ad file"),
    results_path: Path = typer.Option(..., help="Path to NSForest results CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
    organism: str = typer.Option("human", help="Organism: human or mouse"),
):
    """
    Map Ensembl IDs to gene symbols.

    Corresponds to DEMO notebook: Section 9 - Gene mapping
    """
    from commands.gene_mapping import run_gene_mapping
    run_gene_mapping(h5ad_path, results_path, cluster_header, organ, first_author, year, organism)


@app.command("plot-expression")
def plot_expression_command(
    h5ad_path: Path = typer.Option(..., help="Path to input h5ad file"),
    results_path: Path = typer.Option(..., help="Path to NSForest results with gene symbols"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
    plot_types: str = typer.Option("dotplot,violin,matrix", help="Comma-separated: dotplot,violin,matrix"),
    top_n: int = typer.Option(10, help="Number of top markers per cluster"),
):
    """
    Plot expression of top markers (dotplot/violin/matrix).

    Corresponds to DEMO notebook: Section 10 - Expression visualization
    """
    from commands.plot_expression import run_plot_expression
    plots = plot_types.split(',')
    run_plot_expression(h5ad_path, results_path, cluster_header, organ, first_author, year,
                        plots, top_n)


def main():
    app()


if __name__ == "__main__":
    main()
