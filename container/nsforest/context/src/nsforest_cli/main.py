from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence

import scanpy as sc
import typer

# ---- your single-function modules ----
from .dendrogramplot import dendrogramplot_run
from .dotplot import dotplot_run
from .matrixplot import matrixplot_run
from .violinplot import violinplot_run

# existing modules you already have in nsforest_cli/
from .prep_medians import prep_medians_run
from .prep_binary_scores import prep_binary_scores_run
from .run_nsforest import nsforest_run
from .eval_markers import eval_markers_run
from .sanitize import sanitize_labels_run
from .filter_obs import filter_by_obs_run
from .symbolize_genes import symbolize_genes_run
from .build_symbol_map import build_symbol_map_run

app = typer.Typer(add_completion=False)
app = typer.Typer(no_args_is_help=True)

# -----------------------------------------
# Batch lookup the gene symbols -
#  save ENSG as orig - run after all prep,
#  and nsforest but before plotting,
#  optional csv file with mapping
# ----------------------------------------
from .symbolize_genes import symbolize_genes_run
@app.command("build-symbol-map")
def cmd_build_symbol_map(
    *,
    gencode_release: int = typer.Option(..., "--gencode-release", help="Gencode release number (e.g. 49)"),
    out_csv: Path = typer.Option(..., "--out-csv", help="Output CSV path (ensg,symbol)"),
):
    """
    Build ENSG→gene symbol map from Gencode GTF using local curl + unzip.
    """
    from .build_symbol_map import build_symbol_map_run
    build_symbol_map_run(
        gencode_release=gencode_release,
        out_csv=out_csv,
    )
    
@app.command("symbolize")
def cmd_symbolize_genes(
    *,
    h5ad_in: Path = typer.Option(..., "--h5ad-in", help="Input .h5ad file", exists=True, readable=True),
    h5ad_out: Path = typer.Option(..., "--h5ad-out", help="Output .h5ad with gene symbols", dir_okay=False),
    symbol_map_csv: Path = typer.Option(..., "--symbol-map-csv", help="CSV with ENSG,symbol columns", exists=True),
):
    """
    Replace adata.var_names with gene symbols using a local ENSG→symbol CSV.
    Preserves ENSGs in .var['orig_var_names'].
    """
    from .symbolize_genes import symbolize_genes_run
    symbolize_genes_run(
        h5ad_in=h5ad_in,
        h5ad_out=h5ad_out,
        symbol_map_csv=symbol_map_csv,
    )

# --------------------------
# Plotting commands
# --------------------------

# ---------- DENDROGRAM ----------
@app.command("dendrogramplot")
def cmd_dendrogramplot(
    *,
    h5ad_in:     Path                = typer.Option(...,  "--h5ad-in",  help="required h5ad input file ", exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(...,  "--label-key",help="required cluster label"),
    h5ad_out:    Path                = typer.Option(...,  "--h5ad-out", help="required h5ad output file", writable=None, dir_okay=False),
    leaf_range:  Optional[str]       = typer.Option(None, "--leaf-range",   help="Slice of leaf positions, e.g. '0:10'"),
    leaf_indices:Optional[List[int]] = typer.Option(None, "--leaf-indices", help="Explicit leaf indices, e.g. --leaf-indices 0 3 4"),
):
    """
    Display the hierarchical cluster dendrogram.
    ENSG maps to gene symbol if unique.
    Gene symbol names saved to adata.var_names
    nsforest function preprocessing.dendrogram called and saved in format provided.
    """
    adata = dendrogramplot_run(
        h5ad_in=h5ad_in,
        label_key=label_key,
        h5ad_out=h5ad_out,
        leaf_range=leaf_range,
        leaf_indices=leaf_indices,
    )

# ---------- DOTPLOT ----------
@app.command("dotplot")
def cmd_dotplot(
    *,
    h5ad_in:     Path                = typer.Option(..., "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    results_csv: Path                = typer.Option(..., "--results-csv", exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(..., "--label-key", "-l"),
    symbol_map_csv: Optional[Path]   = typer.Option(None, "--symbol-map-csv", help="Optional ENSG→symbol mapping CSV"),
    leaf_range:  Optional[str]       = typer.Option(None, "--leaf-range",   help="Slice of leaf positions, e.g. '0:10'"),
    leaf_indices:Optional[List[int]] = typer.Option(None, "--leaf-indices", help="Explicit leaf indices, e.g. --leaf-indices 0 3 4"),
):
    """
    Display the per-cluster NSForest Markers dotplot (log scaled)
    """
    dotplot_run(
        h5ad_in=h5ad_in,
        results_csv=results_csv,
        label_key=label_key,
        symbol_map_csv=symbol_map_csv,
        leaf_range=leaf_range,
        leaf_indices=leaf_indices,
    )

# ---------- MATRIXPLOT  ----------
@app.command("matrixplot")

def cmd_matrixplot(
    *,
    h5ad_in:     Path                = typer.Option(..., "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    results_csv: Path                = typer.Option(..., "--results-csv", exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(..., "--label-key", "-l"),
    symbol_map_csv: Optional[Path]   = typer.Option(None, "--symbol-map-csv", help="Optional ENSG→symbol mapping CSV"),
    leaf_range:  Optional[str]       = typer.Option(None, "--leaf-range",   help="Slice of leaf positions, e.g. '0:10'"),
    leaf_indices:Optional[List[int]] = typer.Option(None, "--leaf-indices", help="Explicit leaf indices, e.g. --leaf-indices 0 3 4"),
):
    """
    Display the per-cluster NSForest Markers matrixplot (log scaled). Required input nsforest output csv and h5ad file with stored dendrogram.
    """
    matrixplot_run(
        h5ad_in=h5ad_in,
        results_csv=results_csv,
        label_key=label_key,
        symbol_map_csv=symbol_map_csv,
        leaf_range=leaf_range,
        leaf_indices=leaf_indices,
    )

# ---------- VIOLIN ----------
@app.command("violinplot")
def cmd_violinplot(
    *,
    h5ad_in:     Path                = typer.Option(..., "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    results_csv: Path                = typer.Option(..., "--results-csv", exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(..., "--label-key", "-l"),
    symbol_map_csv: Optional[Path]   = typer.Option(None, "--symbol-map-csv", help="Optional ENSG→symbol mapping CSV"),
    leaf_range:  Optional[str]       = typer.Option(None, "--leaf-range",   help="Slice of leaf positions, e.g. '0:10'"),
    leaf_indices:Optional[List[int]] = typer.Option(None, "--leaf-indices", help="Explicit leaf indices, e.g. --leaf-indices 0 3 4"),
):
    """
    Display the per-cluster NSForest Markers violinplot (log scaled). Required input nsforest output csv and h5ad file with stored dendrogram.
    """
    violinplot_run(
        h5ad_in=h5ad_in,
        results_csv=results_csv,
        label_key=label_key,
        symbol_map_csv=symbol_map_csv,
        leaf_range=leaf_range,
        leaf_indices=leaf_indices,
    )

# --------------------------
# NSForest preprocessing / core
# --------------------------
@app.command("prep-medians")
def cmd_prep_medians(
    *,
    h5ad_in:     Path                = typer.Option(...,  "--h5ad-in",  help="required h5ad input file ", exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(...,  "--label-key",help="required cluster label"),
    h5ad_out:    Path                = typer.Option(...,  "--h5ad-out", help="required h5ad output file", writable=None, dir_okay=False),
):
    """
    Compute per-cluster medians (nsforest.pp.prep_medians) and write a new .h5ad.
    """
    adata = prep_medians_run(
                             h5ad_in=h5ad_in,
                             label_key=label_key,
                             h5ad_out=h5ad_out,
                             )


@app.command("prep-binary-scores")
def cmd_prep_binary_scores(
    *,
    h5ad_in:     Path                = typer.Option(...,  "--h5ad-in",  help="required h5ad input file ", exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(...,  "--label-key",help="required cluster label"),
    h5ad_out:    Path                = typer.Option(...,  "--h5ad-out", help="required h5ad output file", writable=None, dir_okay=False),
):
    """
    Compute per-cluster binary scores (nsforest.pp.prep_binary_scores) and write a new .h5ad.
    """
    adata = prep_binary_scores_run(
                                   h5ad_in=h5ad_in,
                                   label_key=label_key,
                                   h5ad_out=h5ad_out,
                                   )


@app.command("nsforest")
def cmd_nsforest(
    *,
    h5ad_in:     Path  = typer.Option(...,  "--h5ad-in",  help="required h5ad input file ", exists=True, dir_okay=False, readable=True),
    label_key:   str   = typer.Option(...,  "--label-key",help="required cluster label"),
    results_csv: Path  = typer.Option(...,  "--results-csv",help="for output with statistics and markers", writable=None, dir_okay=False),
    n_trees:     int   = typer.Option(1000, "--n-trees"),
):
    """
    Run NS-Forest core algorithm and write the results CSV.
    This CLI does not create directories; if output_folder is used, ensure it exists upstream.
    """
    nsforest_run(
        h5ad_in=h5ad_in,
        label_key=label_key,
        results_csv_out=results_csv,
        n_trees=n_trees,
    )


@app.command("eval-markers")
def cmd_eval_markers(
    *,
    h5ad_in:      Path  = typer.Option(...,  "--h5ad-in",    help="required h5ad input file ", exists=True, dir_okay=False, readable=True),
    markers_csv : Path  = typer.Option(...,  "--markers-csv",help="markers per cluster to be evaluated", exists=True, dir_okay=False, readable=True),
    label_key:    str   = typer.Option(...,  "--label-key",  help="required cluster label"),
):
    """
    Evaluate marker sets using nsforest.evaluating helpers; write results CSV.
    """
    eval_markers_run(
        h5ad_in=h5ad_in,
        markers_csv=markers_csv,
        label_key=label_key,
    )


# --------------------------
# Label sanitize & filtering
# --------------------------

@app.command("sanitize-labels")
def cmd_sanitize_labels(
        *,
        h5ad_in:     Path  = typer.Option(...,  "--h5ad-in",  help="required h5ad input file ", exists=True, dir_okay=False, readable=True),
        label_key:   str   = typer.Option(...,  "--label-key",help="required cluster label"),
        h5ad_out:    Path  = typer.Option(...,  "--h5ad-out", help="required h5ad output file", writable=None, dir_okay=False),
):
    """
    Collision-safe label sanitizer (replaces non [A-Za-z0-9_-] with '_', collapses repeats, dedupes).
    Writes a new .h5ad.
    """
    adata = sanitize_labels_run(
            h5ad_in=h5ad_in,
            label_key=label_key,
            h5ad_out=h5ad_out,
    )


@app.command("filter-by-obs")
def cmd_filter_by_obs(
    *,
    h5ad_in:     Path  = typer.Option(...,  "--h5ad-in",  help="required h5ad input file ", exists=True, dir_okay=False, readable=True),
    h5ad_out:    Path  = typer.Option(...,  "--h5ad-out", help="required h5ad output file", writable=None, dir_okay=False),
    obs_key:     str   = typer.Option(..., "--obs-key", "-k", help="obs column to filter"),
    values:      Optional[List[str]] = typer.Option(None, "--values", help="Multiple values (repeat flag)"),
    mode:        str = typer.Option(
        "exact", "--mode", help="exact | contains | regex"
    ),
    case_insensitive: bool = typer.Option(True, "--case-insensitive/--case-sensitive"),
    na_policy: str         = typer.Option("drop", "--na-policy", help="keep | drop"),
    invert: bool           = typer.Option(False, "--invert/--no-invert"),
):
    """
    Filter by a SINGLE obs field. Chain this command multiple times for tissue, then disease, etc.
    Writes the filtered .h5ad to the path you provide.
    """
    if not values:
        typer.echo("[error] You must provide at least one --values entry.", err=True)
        raise typer.Exit(code=1)

    filter_by_obs_run(
        h5ad_in=h5ad_in,
        h5ad_out=h5ad_out,
        obs_key=obs_key,
        values=values,
        mode=mode,
        case_insensitive=case_insensitive,
        na_policy=na_policy,
        invert=invert,
    )


if __name__ == "__main__":
    app()

