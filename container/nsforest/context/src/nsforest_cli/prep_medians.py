# container/nsforest/context/src/nsforest_cli/prep_medians.py
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import typer
from pathlib import Path
import scanpy as sc
import nsforest as ns


def run(
    input_path: Path = typer.Option(..., exists=True, help="Path to input .h5ad file"),
    cluster_header: str = typer.Option("cluster", help="Column name of cluster labels"),
    output_path: Path = typer.Option(..., help="Output path for .h5ad with medians")
):
    """Preprocess AnnData using prep_medians."""
    typer.echo(f"Running prep_medians on: {input_path}")
    adata = sc.read_h5ad(str(input_path))
    adata_subset = ns.pp.prep_medians(adata, cluster_header)
    adata_subset.write_h5ad(str(output_path))
    typer.echo(f"Saved: {output_path}")

