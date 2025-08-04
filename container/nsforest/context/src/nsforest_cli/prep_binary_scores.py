# container/nsforest/context/src/nsforest_cli/prep_binary_scores.py
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import typer
from pathlib import Path
import scanpy as sc
import nsforest as ns


def run(
    input_path: Path = typer.Option(..., exists=True, help="Path to input .h5ad from prep_medians"),
    cluster_header: str = typer.Option("cluster", help="Cluster column header"),
    output_path: Path = typer.Option(..., help="Output path for .h5ad with binary scores")
):
    """Preprocess AnnData using prep_binary_scores."""
    typer.echo(f"Running prep_binary_scores on: {input_path}")
    adata = sc.read_h5ad(str(input_path))
    adata_subset = ns.pp.prep_binary_scores(adata, cluster_header)
    adata_subset.write_h5ad(str(output_path))
    typer.echo(f"Saved: {output_path}")

