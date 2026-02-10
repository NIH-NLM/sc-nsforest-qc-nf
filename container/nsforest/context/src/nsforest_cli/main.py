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
       Corresponds to DEMO notebook: Section 1
    """
    from .dendrogram import run_dendrogram
    run_dendrogram(h5ad_path, cluster_header, organ, first_author, year)

@app.command("cluster-stats")
def cluster_stats_command(
    h5ad_path: Path = typer.Option(..., help="Path to input h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """Compute cluster statistics.
       Corresponds to DEMO notebook: Section 2
    """
    from .cluster_stats import run_cluster_stats
    run_cluster_stats(h5ad_path, cluster_header, organ, first_author, year)

@app.command("prep-medians")
def prep_medians_command(
    h5ad_path: Path = typer.Option(..., help="Path to input h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
    cluster_list: str = typer.Option(None, help="Comma-separated clusters for parallelization"),
):
    """Prepare median expression matrix.
       Corresponds to DEMO notebook: Section 3
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
    Merge partial median CSV files into complete median matrix.
    
    This is the GATHER step after parallelized prep-medians.
    """
    from .merge_medians import run_merge_medians
    files = partial_files.split(',')
    run_merge_medians(files, cluster_header, organ, first_author, year)

@app.command("prep-binary-scores")
def prep_binary_scores_command(
    median_matrix_path: Path = typer.Option(..., help="Path to median matrix CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
    method: str = typer.Option("mwu", help="Binary score method: mwu, ttest"),
):
    """
    Compute binary scores for marker discovery.
    
    Corresponds to DEMO notebook: Section 4 - Prepare binary scores
    Requires complete median matrix from merge-medians step
    """
    from .prep_binary_scores import run_prep_binary_scores
    run_prep_binary_scores(median_matrix_path, cluster_header, organ, first_author, year, method)

def main():
    app()


if __name__ == "__main__":
    main()
