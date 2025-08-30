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

app = typer.Typer(add_completion=False)
app = typer.Typer(no_args_is_help=True)


# --------------------------
# Plotting commands
# --------------------------

# ---------- DOTPLOT ----------
@app.command("dotplot")
def cmd_dotplot(
    h5ad_in:     Path                = typer.Option(..., "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    results_csv: Path                = typer.Option(..., "--results-csv", exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(..., "--label-key", "-l"),
    png_out:     Optional[Path]      = typer.Option(None, "--png-out"),
    svg_out:     Optional[Path]      = typer.Option(None, "--svg-out"),
    leaf_range:  Optional[str]       = typer.Option(None, "--leaf-range",   help="Slice of leaf positions, e.g. '0:10'"),
    leaf_indices:Optional[List[int]] = typer.Option(None, "--leaf-indices", help="Explicit leaf indices, e.g. --leaf-indices 0 3 4"),
):
    """
    Display the per-cluster NSForest Markers dotplot (log scaled)
    """
    if not (png_out or svg_out):
        raise typer.BadParameter("Provide at least one of --png-out / --svg-out")
    dotplot_run(
        h5ad_in=h5ad_in,
        results_csv=results_csv,
        label_key=label_key,
        png_out=png_out,
        svg_out=svg_out,
        leaf_range=leaf_range,
        leaf_indices=leaf_indices,
    )

# ---------- DENDROGRAM ----------
@app.command("dendrogramplot")
def cmd_dendrogramplot(
    h5ad_in:     Path                = typer.Option(...,  "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(...,  "--label-key", "-l"),
    h5ad_out:    Path                = typer.Option(...,  "--h5ad-out",     exists=True, dir_okay=False, readable=True),
    png_out:     Optional[Path]      = typer.Option(None, "--png-out"),
    svg_out:     Optional[Path]      = typer.Option(None, "--svg-out"),
    leaf_range:  Optional[str]       = typer.Option(None, "--leaf-range",   help="Slice of leaf positions, e.g. '0:10'"),
    leaf_indices:Optional[List[int]] = typer.Option(None, "--leaf-indices", help="Explicit leaf indices, e.g. --leaf-indices 0 3 4"),
):
    """
    Display the hierarchical cluster dendrogram - since we are running each step in serial, save the dendrogram output.
    """
    adata = dendrogramplot_run(
        h5ad_in=h5ad_in,
        label_key=label_key,
        png_out=png_out,
        svg_out=svg_out,
        leaf_range=leaf_range,
        leaf_indices=leaf_indices,
    )
    adata.write_h5ad(str(h5ad_out))
   
# ---------- MATRIXPLOT  ----------
@app.command("matrixplot")
def cmd_matrixplot(
    h5ad_in:     Path                = typer.Option(..., "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    results_csv: Path                = typer.Option(..., "--results-csv", exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(..., "--label-key", "-l"),
    png_out:     Optional[Path]      = typer.Option(None, "--png-out"),
    svg_out:     Optional[Path]      = typer.Option(None, "--svg-out"),
    leaf_range:  Optional[str]       = typer.Option(None, "--leaf-range",   help="Slice of leaf positions, e.g. '0:10'"),
    leaf_indices:Optional[List[int]] = typer.Option(None, "--leaf-indices", help="Explicit leaf indices, e.g. --leaf-indices 0 3 4"),
):
    """
    Display the per-cluster NSForest Markers matrixplot (log scaled)
    """
    if not (png_out or svg_out):
        raise typer.BadParameter("Provide at least one of --png-out / --svg-out")
    matrixplot_run(
        h5ad_in=h5ad_in,
        results_csv=results_csv,
        label_key=label_key,
        png_out=png_out,
        svg_out=svg_out,
        leaf_range=leaf_range,
        leaf_indices=leaf_indices,
    )

# ---------- VIOLIN ----------
@app.command("violinplot")
def cmd_violinplot(
    h5ad_in:     Path                = typer.Option(..., "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    results_csv: Path                = typer.Option(..., "--results-csv", exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(..., "--label-key", "-l"),
    png_out:     Optional[Path]      = typer.Option(None, "--png-out"),
    svg_out:     Optional[Path]      = typer.Option(None, "--svg-out"),
    leaf_range:  Optional[str]       = typer.Option(None, "--leaf-range",   help="Slice of leaf positions, e.g. '0:10'"),
    leaf_indices:Optional[List[int]] = typer.Option(None, "--leaf-indices", help="Explicit leaf indices, e.g. --leaf-indices 0 3 4"),
):
    """
    Display the per-cluster NSForest Markers violinplot (log scaled)
    """
    if not (png_out or svg_out):
        raise typer.BadParameter("Provide at least one of --png-out / --svg-out")
    violinplot_run(
        h5ad_in=h5ad_in,
        results_csv=results_csv,
        label_key=label_key,
        png_out=png_out,
        svg_out=svg_out,
        leaf_range=leaf_range,
        leaf_indices=leaf_indices,
    )

# --------------------------
# NSForest preprocessing / core
# --------------------------

@app.command("prep-medians")
def cmd_prep_medians(
    h5ad_in:     Path                = typer.Option(...,  "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(...,  "--label-key", "-l"),
    *,
    h5ad_out:    Path                = typer.Option(...,  "--h5ad-out",     exists=True, dir_okay=False, readable=True),
):
    """
    Compute per-cluster medians (nsforest.pp.prep_medians) and write a new .h5ad.
    """
    adata = prep_medians_run(h5ad_in=h5ad_in,
                             label_key=label_key)
    adata.write_h5ad(str(h5ad_out))


@app.command("prep-binary-scores")
def cmd_prep_binary_scores(
    h5ad_in:     Path                = typer.Option(...,  "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(...,  "--label-key", "-l"),
    *,
    h5ad_out:    Path                = typer.Option(...,  "--h5ad-out",     exists=True, dir_okay=False, readable=True),
):
    """
    Compute per-cluster binary scores (nsforest.pp.prep_binary_scores) and write a new .h5ad.
    """
    adata = prep_binary_scores_run(h5ad_in=h5ad_in,
                                   label_key=label_key)
    adata.write_h5ad(str(h5ad_out))


@app.command("nsforest")
def cmd_nsforest(
    h5ad_in:     Path                = typer.Option(...,  "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    results_csv: Path                = typer.Option(...,  "--results-csv", exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(...,  "--label-key", "-l"),
    h5ad_out:    Path                = typer.Option(...,  "--h5ad-out",     exists=True, dir_okay=False, readable=True),
    *,
    n_trees: int                     = typer.Option(1000, "--n-trees"),
):
    """
    Run NS-Forest core algorithm and write the results CSV.
    This CLI does not create directories; if output_folder is used, ensure it exists upstream.
    """
    run_nsforest_run(
        h5ad_in=h5ad_in,
        label_key=label_key,
        results_csv_out=results_csv_out,
        n_trees=n_trees,
    )


@app.command("eval-markers")
def cmd_eval_markers(
    h5ad_in:     Path                = typer.Option(...,  "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    eval_csv_out: Path               = typer.Option(...,  "--eval-results-csv", exists=True, dir_okay=False, readable=True),
    *,
    label_key:   str                 = typer.Option(...,  "--label-key", "-l"),
):
    """
    Evaluate marker sets using nsforest.evaluating helpers; write results CSV.
    """
    eval_markers_run(
        h5ad_in=h5ad_in,
        eval_csv_out=eval_csv_out,
        label_key=label_key
    )


# --------------------------
# Label sanitize & filtering
# --------------------------

@app.command("sanitize-labels")
def cmd_sanitize_labels(
    h5ad_in:     Path                = typer.Option(...,  "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(...,  "--label-key", "-l"),
    *,
    h5ad_out:    Path                = typer.Option(...,  "--h5ad-out",     exists=True, dir_okay=False, readable=True),
):
    """
    Collision-safe label sanitizer (replaces non [A-Za-z0-9_-] with '_', collapses repeats, dedupes).
    Writes a new .h5ad to the path you provide.
    """
    adata = sanitize_labels_run(
        h5ad_in=h5ad_in,
        label_key=label_key)
    adata.write_h5ad(str(h5ad_out))


@app.command("filter-by-obs")
def cmd_filter_by_obs(
    h5ad_in:     Path                = typer.Option(...,  "--h5ad-in",     exists=True, dir_okay=False, readable=True),
    label_key:   str                 = typer.Option(...,  "--label-key", "-l"),
    h5ad_out:    Path                = typer.Option(...,  "--h5ad-out",     exists=True, dir_okay=False, readable=True),
    obs_key:     str                 = typer.Option(..., "--obs-key", "-k", help="obs column to filter"),
    value:       Optional[str]       = typer.Option(None, "--value", "-v", help="Single exact value"),
    values:      Optional[List[str]] = typer.Option(None, "--values", help="Multiple values (repeat flag)"),
    mode:        str = typer.Option(
        "exact", "--mode", help="exact | contains | regex"
    ),
    case_insensitive: bool           = typer.Option(True, "--case-insensitive/--case-sensitive"),
    na_policy: str                   = typer.Option("drop", "--na-policy", help="keep | drop"),
    *,
    invert: bool                     = typer.Option(False, "--invert/--no-invert"),
):
    """
    Filter by a SINGLE obs field. Chain this command multiple times for tissue, then disease, etc.
    Writes the filtered .h5ad to the path you provide.
    """
    adata = filter_by_obs_run(
        h5ad_in=h5ad_in,
        obs_key=obs_key,
        values=[value] if (value and not values) else values,
        mode=mode,
        case_insensitive=case_insensitive,
        na_policy=na_policy,
        invert=invert,
    )
    adata.write_h5ad(str(h5ad_out))


if __name__ == "__main__":
    app()

