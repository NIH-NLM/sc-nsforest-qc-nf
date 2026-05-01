"""
NSForest CLI - Command-line interface for NSForest workflow.

Modular commands matching DEMO_NS-Forest_workflow.py workflow.
"""

import typer
from pathlib import Path

app = typer.Typer(help="NSForest workflow commands")

@app.command("cluster-cid-mapping")
def cluster_cid_mapping_command(
    h5ad_path: Path = typer.Option(..., help="Path to h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    journal: str = typer.Option(..., help="Journal"),
    year: str = typer.Option(..., help="Publication year"),
    embedding: str = typer.Option("", help="Embedding key"),
    dataset_version_id: str = typer.Option("", help="Dataset version ID"),
    cid_column: str = typer.Option(
        "cell_type_ontology_term_id",
        help="obs column holding current cell-ontology term ID",
    ),
    default_skos: str = typer.Option(
        "",
        help='Value to pre-populate in skos column: "", "exact", or "related"',
    ),
):
    """Emit cluster -> cell_ontology_id mapping CSV (cluster_name, skos, manual_mapped_cid, cell_ontology_id)."""
    from .cluster_cid_mapping import run_cluster_cid_mapping
    run_cluster_cid_mapping(
        h5ad_path, cluster_header, organ, first_author, journal, year,
        embedding, dataset_version_id,
        cid_column=cid_column,
        default_skos=default_skos,
    )

@app.command("cluster-stats")
def cluster_stats_command(
    h5ad_path: Path = typer.Option(..., help="Path to h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    journal: str = typer.Option(..., help="Journal"),
    year: str = typer.Option(..., help="Publication year"),
    embedding: str = typer.Option("", help="Embedding key"),
    dataset_version_id: str = typer.Option("", help="Dataset version ID"),

):
    """Compute cluster statistics."""
    from .cluster_stats import run_cluster_stats
    run_cluster_stats(h5ad_path, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id)

@app.command("dendrogram")
def dendrogram_command(
    h5ad_path: Path = typer.Option(..., help="Path to h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    journal: str = typer.Option(..., help="Journal"),
    year: str = typer.Option(..., help="Publication year"),
    embedding: str = typer.Option("", help="Embedding key"),
    dataset_version_id: str = typer.Option("", help="Dataset version ID"),

):
    """Generate dendrogram and cluster statistics."""
    from .dendrogram import run_dendrogram
    run_dendrogram(h5ad_path, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id)

@app.command("filter-adata")
def filter_adata_command(
    h5ad_path: Path = typer.Option(..., help="Path to h5ad file"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    journal: str = typer.Option(..., help="Journal"),
    year: str = typer.Option(..., help="Publication year"),
    embedding: str = typer.Option("", help="Embedding key"),
    dataset_version_id: str = typer.Option("", help="Dataset version ID"),
    filter_normal: bool = typer.Option(True, help="Filter to normal adult cells only"),
    uberon: Path = typer.Option(None, help="UBERON JSON from cellxgene-harvester resolve-uberon"),
    disease: Path = typer.Option(None, help="Disease JSON from cellxgene-harvester resolve-disease"),
    hsapdv: Path = typer.Option(None, help="HsapDv JSON from cellxgene-harvester resolve-hsapdv --min-age N"),
    min_cluster_size: int = typer.Option(5, help="Minimum cells per cluster"),
    tissue_ontology_term_id: str = typer.Option(None, help="Pipe-separated UBERON term IDs from CSV row"),
    disease_ontology_term_id: str = typer.Option(None, help="Pipe-separated disease ontology term IDs from CSV row"),
    development_stage_ontology_term_id: str = typer.Option(None, help="Pipe-separated HsapDv term IDs from CSV row"),
    filter_obs_column: str = typer.Option('', help="obs column to filter on (e.g., Dataset)"),
    filter_obs_value: str = typer.Option('', help="Value to keep in filter_obs_column (e.g., 'Tosti et al. 2021')"),

):
    """Filter adata by tissue, disease, age, and minimum cluster size."""
    from .filter_adata import run_filter_adata
    run_filter_adata(
        h5ad_path          = h5ad_path,
        cluster_header     = cluster_header,
        organ              = organ,
        first_author       = first_author,
        journal            = journal,
        year               = year,
        embedding          = embedding,
        dataset_version_id = dataset_version_id,
        filter_normal      = filter_normal,
        filter_obs_column  = filter_obs_column or None,
        filter_obs_value   = filter_obs_value or None,
        uberon_json        = str(uberon)  if uberon  else None,
        disease_json       = str(disease) if disease else None,
        hsapdv_json        = str(hsapdv)  if hsapdv  else None,
        min_cluster_size   = min_cluster_size,
        row_uberon_ids     = tissue_ontology_term_id,
        row_disease_ids    = disease_ontology_term_id,
        row_hsapdv_ids     = development_stage_ontology_term_id,
    )

@app.command("generate-s3-manifest")
def generate_s3_manifest_command(
    s3_base: str = typer.Option(..., help="S3 results base path (e.g. s3://bucket/.../jobs/{id}/results)"),
):
    """Generate master_s3_manifest.csv listing all output files and their S3 paths."""
    from .generate_s3_manifest import run_generate_s3_manifest
    run_generate_s3_manifest(s3_base)

@app.command("merge-nsforest-results")
def merge_nsforest_results_command(
    partial_files: str = typer.Option(..., help="Comma-separated list of partial results CSV files"),
    filtered_h5ad: Path = typer.Option(..., help="Path to adata_filtered.h5ad (for gene_symbol mapping)"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    journal: str = typer.Option(..., help="Journal"),
    year: str = typer.Option(..., help="Publication year"),
    embedding: str = typer.Option("", help="Embedding key"),
    dataset_version_id: str = typer.Option("", help="Dataset version ID"),

):
    """Merge partial NSForest results files and save csv + pkl."""
    from .merge_nsforest_results import run_merge_nsforest_results
    files = partial_files.split(',')
    run_merge_nsforest_results(files, filtered_h5ad, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id)

@app.command("plot-histograms")
def plot_histograms_command(
    medians_csv: Path = typer.Option(..., help="Path to medians CSV"),
    binary_scores_csv: Path = typer.Option(..., help="Path to binary scores CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    journal: str = typer.Option(..., help="Journal"),
    year: str = typer.Option(..., help="Publication year"),
    embedding: str = typer.Option("", help="Embedding key"),
    dataset_version_id: str = typer.Option("", help="Dataset version ID"),

):
    """Plot histograms of non-zero median and binary score values."""
    from .plot_histograms import run_plot_histograms
    run_plot_histograms(medians_csv, binary_scores_csv, cluster_header, organ, first_author, journal, year,
                        embedding, dataset_version_id)

@app.command("plots")
def plots_command(
    h5ad_path: Path = typer.Option(..., help="Path to adata_filtered.h5ad"),
    results_csv: Path = typer.Option(..., help="Path to NSForest results CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    journal: str = typer.Option(..., help="Journal"),
    year: str = typer.Option(..., help="Publication year"),
    embedding: str = typer.Option("", help="Embedding key"),
    dataset_version_id: str = typer.Option("", help="Dataset version ID"),

):
    """Create NSForest visualization plots with gene symbol mapping."""
    from .plots import run_plots
    run_plots(h5ad_path, results_csv, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id)

@app.command("prep-binary-scores")
def prep_binary_scores_command(
    h5ad_path: Path = typer.Option(..., help="Path to adata_filtered.h5ad"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    journal: str = typer.Option(..., help="Journal"),
    year: str = typer.Option(..., help="Publication year"),
    embedding: str = typer.Option("", help="Embedding key"),
    dataset_version_id: str = typer.Option("", help="Dataset version ID"),

):
    """Compute binary scores per cluster. Saves binary_scores csv + pkl."""
    from .prep_binary_scores import run_prep_binary_scores
    run_prep_binary_scores(h5ad_path, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id)


@app.command("prep-medians")
def prep_medians_command(
    h5ad_path: Path = typer.Option(..., help="Path to adata_filtered.h5ad"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    journal: str = typer.Option(..., help="Journal"),
    year: str = typer.Option(..., help="Publication year"),
    embedding: str = typer.Option("", help="Embedding key"),
    dataset_version_id: str = typer.Option("", help="Dataset version ID"),

):
    """Compute median expression per cluster. Saves medians csv + pkl."""
    from .prep_medians import run_prep_medians
    run_prep_medians(h5ad_path, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id)

@app.command("run-nsforest")
def run_nsforest_command(
    h5ad_path: Path = typer.Option(..., help="Path to adata_filtered.h5ad"),
    medians_csv: Path = typer.Option(..., help="Path to medians CSV"),
    binary_scores_csv: Path = typer.Option(..., help="Path to binary scores CSV"),
    cluster_header: str = typer.Option(..., help="Column name for clusters"),
    organ: str = typer.Option(..., help="Organ/tissue"),
    first_author: str = typer.Option(..., help="First author"),
    journal: str = typer.Option(..., help="Journal"),
    year: str = typer.Option(..., help="Publication year"),
    embedding: str = typer.Option("", help="Embedding key"),
    dataset_version_id: str = typer.Option("", help="Dataset version ID"),
    cluster_list: str = typer.Option(None, help="Comma-separated cluster list (for parallelization)"),
    n_trees: int = typer.Option(1000, help="Number of trees in random forest"),
    n_genes_eval: int = typer.Option(6, help="Number of top genes to evaluate"),
):
    """Run NSForest algorithm to identify marker genes."""
    from .run_nsforest import run_nsforest
    clusters = cluster_list.split(',') if cluster_list else None
    run_nsforest(h5ad_path, medians_csv, binary_scores_csv, cluster_header,
                 organ, first_author, journal, year, embedding, dataset_version_id, clusters, n_trees, n_genes_eval)


if __name__ == "__main__":
    app()
