# repo: container/nsforest/context/src/nsforest_cli
# change set: add dendrogram writer (png+svg+h5ad), fix dotplot to log-scale, add matrixplot; ensure dendrogram first; optional gene symbol conversion via REST util

# =========================
# FILE: nsforest_cli/plotting_ext.py
# =========================
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence, Union, Optional, Callable

import scanpy as sc  # type: ignore
from anndata import AnnData  # type: ignore

# NS-Forest public API
# - preprocessing lives under ns.pp.*
# - plotting lives under ns.pl.*
#   (refs: https://nsforest.readthedocs.io/en/latest/preprocessing.html#preprocessing.dendrogram
#          https://nsforest.readthedocs.io/en/latest/plotting_scanpy.html)
from nsforest import ns, utils as ns_utils  # type: ignore


# ---- types ----
Markers = Union[Sequence[str], Dict[str, Sequence[str]]]
GeneResolver = Callable[[Sequence[str], str], Sequence[str]]


def _ensure_outdir(path: Union[str, Path]) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def resolve_gene_symbols(
    genes: Sequence[str],
    species: str = "human",
    resolver: Optional[GeneResolver] = None,
) -> List[str]:
    """Resolve input gene identifiers to gene symbols via existing REST utility.

    Why: keep plots consistent with downstream steps that expect gene symbols.
    Falls back to identity if resolver is unavailable.
    """
    if resolver is not None:
        try:
            return list(resolver(list(genes), species))
        except Exception:
            # keep moving if remote service hiccups
            pass

    # Try dynamic import of an existing repo utility if present
    try:
        # expected signature: to_gene_symbols_via_rest(genes: List[str], species: str) -> List[str]
        from nsforest_cli.utils import to_gene_symbols_via_rest  # type: ignore

        return list(to_gene_symbols_via_rest(list(genes), species))
    except Exception:
        # No-op fallback
        return list(genes)


def save_dendrogram_and_h5ad(
    adata: AnnData,
    cluster_header: str,
    *,
    output_folder: Union[str, Path],
    outputfilename_suffix: str = "",
    h5ad_out: Optional[Union[str, Path]] = None,
    **kwargs,
) -> Path:
    """Generate scanpy/NS-Forest dendrogram and persist to disk.

    - Adds `adata.uns["dendrogram_{cluster_header}"]` in-memory
    - Saves PNG + SVG dendrogram images
    - Writes a new h5ad containing the dendrogram in `.uns`
    """
    outdir = _ensure_outdir(output_folder)

    # NS-Forest wrapper delegates to scanpy.tl/pl.dendrogram
    # (docs: preprocessing.dendrogram(..., save='png' | 'svg'))
    ns.pp.dendrogram(
        adata,
        cluster_header,
        plot=True,
        save="png",
        output_folder=str(outdir),
        outputfilename_suffix=outputfilename_suffix,
        **kwargs,
    )
    ns.pp.dendrogram(
        adata,
        cluster_header,
        plot=True,
        save="svg",
        output_folder=str(outdir),
        outputfilename_suffix=outputfilename_suffix,
        **kwargs,
    )

    # persist h5ad with dendrogram in .uns
    if h5ad_out is None:
        h5ad_out = outdir / f"{cluster_header}_with_dendrogram.h5ad"
    h5ad_out = Path(h5ad_out)
    h5ad_out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(h5ad_out)
    return h5ad_out


def scaled_dotplot(
    adata: AnnData,
    markers: Markers,
    cluster_header: str,
    *,
    output_folder: Union[str, Path],
    outputfilename_suffix: str = "",
    dendrogram: Union[bool, List[str]] = True,
    log: bool = True,
    **kwargs,
) -> None:
    """Save PNG + SVG dotplots (log-scaled) using NS-Forest wrapper.

    Important: `dendrogram=True` tells NS-Forest/scanpy to use
    `adata.uns["dendrogram_{cluster_header}"]` ordering.
    """
    outdir = _ensure_outdir(output_folder)

    ns.pl.dotplot(
        adata,
        markers,
        cluster_header,
        dendrogram=dendrogram,
        log=log,
        save="png",
        output_folder=str(outdir),
        outputfilename_suffix=outputfilename_suffix,
        **kwargs,
    )
    ns.pl.dotplot(
        adata,
        markers,
        cluster_header,
        dendrogram=dendrogram,
        log=log,
        save="svg",
        output_folder=str(outdir),
        outputfilename_suffix=outputfilename_suffix,
        **kwargs,
    )


def matrixplot(
    adata: AnnData,
    markers: Markers,
    cluster_header: str,
    *,
    output_folder: Union[str, Path],
    outputfilename_suffix: str = "",
    dendrogram: Union[bool, List[str]] = True,
    log: bool = True,
    **kwargs,
) -> None:
    """Save PNG + SVG matrixplots (log-scaled) using NS-Forest wrapper."""
    outdir = _ensure_outdir(output_folder)

    ns.pl.matrixplot(
        adata,
        markers,
        cluster_header,
        dendrogram=dendrogram,
        log=log,
        save="png",
        output_folder=str(outdir),
        outputfilename_suffix=outputfilename_suffix,
        **kwargs,
    )
    ns.pl.matrixplot(
        adata,
        markers,
        cluster_header,
        dendrogram=dendrogram,
        log=log,
        save="svg",
        output_folder=str(outdir),
        outputfilename_suffix=outputfilename_suffix,
        **kwargs,
    )


# =========================
# FILE: nsforest_cli/main.py  (patched)
# =========================
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import pandas as pd  # type: ignore
import scanpy as sc  # type: ignore
import typer

from nsforest import utils as ns_utils  # type: ignore

from .plotting_ext import (
    save_dendrogram_and_h5ad,
    scaled_dotplot,
    matrixplot,
    resolve_gene_symbols,
)

app = typer.Typer(help="NSForest CLI: dendrogram + plots orchestrator")


@app.command("prep-dendrogram")
def prep_dendrogram(
    h5ad: Path = typer.Argument(..., exists=True, readable=True, help="Input .h5ad"),
    cluster_header: str = typer.Option(..., "--cluster-header", "-c", help="obs column for clusters"),
    output_folder: Path = typer.Option(Path("outputs"), "--out", "-o"),
    outputfilename_suffix: str = typer.Option("", "--suffix"),
    h5ad_out: Optional[Path] = typer.Option(None, "--h5ad-out", help="Write new .h5ad here"),
):
    """Generate dendrogram first, persist to disk (PNG+SVG+h5ad)."""
    adata = sc.read_h5ad(str(h5ad))
    typer.echo(f"Loaded AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

    out_h5ad = save_dendrogram_and_h5ad(
        adata,
        cluster_header,
        output_folder=output_folder,
        outputfilename_suffix=outputfilename_suffix or cluster_header,
        h5ad_out=h5ad_out,
    )
    typer.echo(f"Saved: {out_h5ad}")


@app.command("plot-markers")
def plot_markers(
    h5ad: Path = typer.Argument(..., exists=True, readable=True, help=".h5ad with dendrogram in .uns"),
    cluster_header: str = typer.Option(..., "--cluster-header", "-c"),
    markers_csv: Optional[Path] = typer.Option(None, "--markers-csv", help="CSV with columns [clusterName, markers]"),
    col_cluster: str = typer.Option("clusterName", help="CSV column: cluster"),
    col_markers: str = typer.Option("markers", help="CSV column: markers (list or pipe/csv)"),
    species: str = typer.Option("human", help="Species for gene symbol resolver"),
    output_folder: Path = typer.Option(Path("outputs"), "--out", "-o"),
    outputfilename_suffix: str = typer.Option("", "--suffix"),
    save_violin: bool = typer.Option(True, help="Also save stacked violin plots"),
):
    """Plot dotplot (log), matrixplot (log), and stackedviolin using dendrogram order.

    Order enforced: dendrogram must exist in `adata.uns` before plotting.
    """
    adata = sc.read_h5ad(str(h5ad))

    # Ensure dendrogram exists
    dendro_key = f"dendrogram_{cluster_header}"
    if dendro_key not in adata.uns:
        # generate quickly without saving again (keeps order consistent)
        ns.pp.dendrogram(adata, cluster_header, plot=False)

    # Prepare markers dict
    if markers_csv is None:
        raise typer.BadParameter("--markers-csv is required")

    df = pd.read_csv(markers_csv)
    markers_dict = ns_utils.prepare_markers(
        df, col_cluster, col_markers, output_folder=str(output_folder), outputfilename_prefix=outputfilename_suffix or cluster_header
    )

    # Convert each gene list to symbols via REST utility (if available)
    resolved: Dict[str, List[str]] = {}
    for clust, genes in markers_dict.items():
        resolved[clust] = resolve_gene_symbols(genes, species=species)

    # Save dotplot (log) PNG + SVG
    scaled_dotplot(
        adata,
        resolved,
        cluster_header,
        output_folder=output_folder,
        outputfilename_suffix=outputfilename_suffix or cluster_header,
        dendrogram=True,
        log=True,
    )

    # Save violin PNG + SVG (NS-Forest stacked violin)
    if save_violin:
        # NS-Forest wrapper exposes save=<fmt>
        ns.pl.stackedviolin(
            adata,
            resolved,
            cluster_header,
            dendrogram=True,
            save="png",
            output_folder=str(output_folder),
            outputfilename_suffix=outputfilename_suffix or cluster_header,
        )
        ns.pl.stackedviolin(
            adata,
            resolved,
            cluster_header,
            dendrogram=True,
            save="svg",
            output_folder=str(output_folder),
            outputfilename_suffix=outputfilename_suffix or cluster_header,
        )

    # Save matrixplot (log) PNG + SVG
    matrixplot(
        adata,
        resolved,
        cluster_header,
        output_folder=output_folder,
        outputfilename_suffix=outputfilename_suffix or cluster_header,
        dendrogram=True,
        log=True,
    )


if __name__ == "__main__":
    app()

