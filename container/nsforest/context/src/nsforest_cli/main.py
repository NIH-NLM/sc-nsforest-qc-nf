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


def main():
    app()


if __name__ == "__main__":
    main()
