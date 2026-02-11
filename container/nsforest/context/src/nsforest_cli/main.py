"""
NSForest CLI - Command-line interface for NSForest workflow.

Modular commands matching DEMO_NS-forest_workflow.py workflow.
"""

import typer
from pathlib import Path

app = typer.Typer(help="NSForest workflow commands")


@app.command("filter-adata")
def filter_adata_command(
    h5ad_path: Path = typer.Option(..., help="Path to h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
    filter_normal: bool = typer.Option(False, help="Filter to normal cells only"),
    tissue: str = typer.Option(None, help="Target tissue to filter"),
    disease_column: str = typer.Option("disease", help="Disease column name in adata.obs"),
    tissue_column: str = typer.Option("tissue", help="Tissue column name in adata.obs"),
    min_cluster_size: int = typer.Option(5, help="Minimum cells per cluster"),
):
    """
    Filter adata by disease, tissue, and minimum cluster size.
    
    Creates before/after dendrograms and statistics.
    Outputs filtered h5ad file.
    """
    from .filter_adata import run_filter_adata
    run_filter_adata(h5ad_path, cluster_header, organ, first_author, year,
                     filter_normal, tissue, disease_column, tissue_column, min_cluster_size)


@app.command("dendrogram")
def dendrogram_command(
    h5ad_path: Path = typer.Option(..., help="Path to h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """
    Generate dendrogram and cluster statistics.
    
    Creates: dendrogram.svg, cluster_sizes.csv, cluster_order.csv, summary_normal.csv
    """
    from .dendrogram import run_dendrogram
    run_dendrogram(h5ad_path, cluster_header, organ, first_author, year)


@app.command("prep-medians")
def prep_medians_command(
    h5ad_path: Path = typer.Option(..., help="Path to h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
    cluster_list: str = typer.Option(None, help="Comma-separated cluster list (for parallelization)"),
):
    """
    Compute median expression per cluster.
    
    Uses ns.pp.prep_medians() to filter positive genes and compute medians.
    """
    from .prep_medians import run_prep_medians
    clusters = cluster_list.split(',') if cluster_list else None
    run_prep_medians(h5ad_path, cluster_header, organ, first_author, year, clusters)


@app.command("merge-medians")
def merge_medians_command(
    partial_files: str = typer.Option(..., help="Comma-separated list of partial median CSV files"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """
    Merge partial median files and save csv + pkl.
    """
    from .merge_medians import run_merge_medians
    files = partial_files.split(',')
    run_merge_medians(files, cluster_header, organ, first_author, year)


@app.command("prep-binary-scores")
def prep_binary_scores_command(
    h5ad_path: Path = typer.Option(..., help="Path to h5ad file (adata_prep from prep_medians)"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
    cluster_list: str = typer.Option(None, help="Comma-separated cluster list (for parallelization)"),
):
    """
    Compute binary scores per cluster.
    
    Uses ns.pp.prep_binary_scores() on filtered adata from prep_medians.
    """
    from .prep_binary_scores import run_prep_binary_scores
    clusters = cluster_list.split(',') if cluster_list else None
    run_prep_binary_scores(h5ad_path, cluster_header, organ, first_author, year, clusters)


@app.command("merge-binary-scores")
def merge_binary_scores_command(
    partial_files: str = typer.Option(..., help="Comma-separated list of partial binary scores CSV files"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """
    Merge partial binary scores files and save csv + pkl.
    """
    from .merge_binary_scores import run_merge_binary_scores
    files = partial_files.split(',')
    run_merge_binary_scores(files, cluster_header, organ, first_author, year)


@app.command("plot-histograms")
def plot_histograms_command(
    medians_csv: Path = typer.Option(..., help="Path to medians CSV"),
    binary_scores_csv: Path = typer.Option(..., help="Path to binary scores CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """
    Plot histograms of non-zero median and binary score values.
    """
    from .plot_histograms import run_plot_histograms
    run_plot_histograms(medians_csv, binary_scores_csv, cluster_header, organ, first_author, year)


@app.command("run-nsforest")
def run_nsforest_command(
    h5ad_path: Path = typer.Option(..., help="Path to adata_prep h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
    cluster_list: str = typer.Option(None, help="Comma-separated cluster list (for parallelization)"),
    n_trees: int = typer.Option(1000, help="Number of trees in random forest"),
    n_genes_eval: int = typer.Option(6, help="Number of top genes to evaluate"),
):
    """
    Run NSForest algorithm to identify marker genes.
    
    Uses nsforesting.NSForest() with cluster_list for parallelization.
    """
    from .run_nsforest import run_nsforest
    clusters = cluster_list.split(',') if cluster_list else None
    run_nsforest(h5ad_path, cluster_header, organ, first_author, year, clusters, n_trees, n_genes_eval)


@app.command("merge-nsforest-results")
def merge_nsforest_results_command(
    partial_files: str = typer.Option(..., help="Comma-separated list of partial results CSV files"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """
    Merge partial NSForest results files and save csv + pkl.
    """
    from .merge_nsforest_results import run_merge_nsforest_results
    files = partial_files.split(',')
    run_merge_nsforest_results(files, cluster_header, organ, first_author, year)


@app.command("plots")
def plots_command(
    h5ad_path: Path = typer.Option(..., help="Path to original h5ad file"),
    results_csv: Path = typer.Option(..., help="Path to NSForest results CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """
    Create NSForest visualization plots with gene symbol mapping.
    
    Generates boxplots, scatter plots, and expression plots (dotplot, violin, matrix).
    """
    from .plots import run_plots
    run_plots(h5ad_path, results_csv, cluster_header, organ, first_author, year)


if __name__ == "__main__":
    app()
