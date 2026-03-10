"""
NSForest CLI - Command-line interface for NSForest workflow.

Modular commands matching DEMO_NS-Forest_workflow.py workflow.
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
    filter_normal: bool = typer.Option(True, help="Filter to normal adult cells only"),
    uberon: Path = typer.Option(None, help="UBERON JSON from cellxgene-harvester resolve-uberon"),
    disease: Path = typer.Option(None, help="Disease JSON from cellxgene-harvester resolve-disease"),
    hsapdv: Path = typer.Option(None, help="HsapDv JSON from cellxgene-harvester resolve-hsapdv --min-age N"),
    min_cluster_size: int = typer.Option(5, help="Minimum cells per cluster"),
    tissue_ontology_term_id: str = typer.Option(None, help="Pipe-separated UBERON term IDs from CSV row"),
    disease_ontology_term_id: str = typer.Option(None, help="Pipe-separated disease ontology term IDs from CSV row"),
    development_stage_ontology_term_id: str = typer.Option(None, help="Pipe-separated HsapDv term IDs from CSV row"),
):
    """Filter adata by tissue, disease, age, and minimum cluster size."""
    from .filter_adata import run_filter_adata
    run_filter_adata(
        h5ad_path        = h5ad_path,
        cluster_header   = cluster_header,
        organ            = organ,
        first_author     = first_author,
        year             = year,
        filter_normal    = filter_normal,
        uberon_json      = str(uberon)  if uberon  else None,
        disease_json     = str(disease) if disease else None,
        hsapdv_json      = str(hsapdv)  if hsapdv  else None,
        min_cluster_size = min_cluster_size,
        row_uberon_ids   = tissue_ontology_term_id,
        row_disease_ids  = disease_ontology_term_id,
        row_hsapdv_ids   = development_stage_ontology_term_id,
    )


@app.command("dendrogram")
def dendrogram_command(
    h5ad_path: Path = typer.Option(..., help="Path to h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """Generate dendrogram and cluster statistics."""
    from .dendrogram import run_dendrogram
    run_dendrogram(h5ad_path, cluster_header, organ, first_author, year)


@app.command("prep-medians")
def prep_medians_command(
    h5ad_path: Path = typer.Option(..., help="Path to adata_filtered.h5ad"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """Compute median expression per cluster. Saves medians csv + pkl."""
    from .prep_medians import run_prep_medians
    run_prep_medians(h5ad_path, cluster_header, organ, first_author, year)


@app.command("prep-binary-scores")
def prep_binary_scores_command(
    h5ad_path: Path = typer.Option(..., help="Path to adata_filtered.h5ad"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """Compute binary scores per cluster. Saves binary_scores csv + pkl."""
    from .prep_binary_scores import run_prep_binary_scores
    run_prep_binary_scores(h5ad_path, cluster_header, organ, first_author, year)


@app.command("plot-histograms")
def plot_histograms_command(
    medians_csv: Path = typer.Option(..., help="Path to medians CSV"),
    binary_scores_csv: Path = typer.Option(..., help="Path to binary scores CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """Plot histograms of non-zero median and binary score values."""
    from .plot_histograms import run_plot_histograms
    run_plot_histograms(medians_csv, binary_scores_csv, cluster_header, organ, first_author, year)


@app.command("run-nsforest")
def run_nsforest_command(
    h5ad_path: Path = typer.Option(..., help="Path to adata_filtered.h5ad"),
    medians_csv: Path = typer.Option(..., help="Path to medians CSV"),
    binary_scores_csv: Path = typer.Option(..., help="Path to binary scores CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
    cluster_list: str = typer.Option(None, help="Comma-separated cluster list (for parallelization)"),
    n_trees: int = typer.Option(1000, help="Number of trees in random forest"),
    n_genes_eval: int = typer.Option(6, help="Number of top genes to evaluate"),
):
    """Run NSForest algorithm to identify marker genes."""
    from .run_nsforest import run_nsforest
    clusters = cluster_list.split(',') if cluster_list else None
    run_nsforest(h5ad_path, medians_csv, binary_scores_csv, cluster_header,
                 organ, first_author, year, clusters, n_trees, n_genes_eval)


@app.command("merge-nsforest-results")
def merge_nsforest_results_command(
    partial_files: str = typer.Option(..., help="Comma-separated list of partial results CSV files"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """Merge partial NSForest results files and save csv + pkl."""
    from .merge_nsforest_results import run_merge_nsforest_results
    files = partial_files.split(',')
    run_merge_nsforest_results(files, cluster_header, organ, first_author, year)


@app.command("plots")
def plots_command(
    h5ad_path: Path = typer.Option(..., help="Path to adata_filtered.h5ad"),
    results_csv: Path = typer.Option(..., help="Path to NSForest results CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    year: str = typer.Option(..., help="Publication year"),
):
    """Create NSForest visualization plots with gene symbol mapping."""
    from .plots import run_plots
    run_plots(h5ad_path, results_csv, cluster_header, organ, first_author, year)


if __name__ == "__main__":
    app()
