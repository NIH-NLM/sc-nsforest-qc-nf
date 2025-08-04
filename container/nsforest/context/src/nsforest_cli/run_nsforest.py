# container/nsforest/context/src/nsforest_cli/run_nsforest.py
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import typer
from pathlib import Path
import scanpy as sc
import nsforesting


def run(
    input_path: Path = typer.Option(..., exists=True, help="Path to input .h5ad after binary scores"),
    cluster_header: str = typer.Option("cluster", help="Cluster column header"),
    output_folder: Path = typer.Option(Path("./NSForest_outputs"), help="Output directory")
):
    """Run NS-Forest marker gene selection."""
    typer.echo(f"Running NSForest on: {input_path}")
    adata = sc.read_h5ad(str(input_path))
    output_folder.mkdir(parents=True, exist_ok=True)

    nsforesting.NSForest(
        adata,
        cluster_header,
        save=True,
        save_supplementary=True,
        output_folder=str(output_folder),
        outputfilename_prefix=cluster_header
    )

    typer.echo(f"NSForest complete. Output saved to: {output_folder}")

