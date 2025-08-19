# container/nsforest/context/src/nsforest_cli/main.py
from __future__ import annotations

from pathlib import Path
from typing import List, Optional
import typer
import plotly.io as pio

# ---- pure helpers: each provides exactly one function ----
from .dotplot            import dotplot_run
from .violinplot         import violinplot_run
from .dendrogramplot     import dendrogramplot_run
from .dendro_subset      import subset_by_dendro_range_run, print_dendro_order_run

from .prep_medians       import prep_medians_run
from .prep_binary_scores import prep_binary_scores_run
from .run_nsforest       import nsforest_run
from .eval_markers       import eval_markers_run

from .sanitize           import sanitize_labels_run
from .filter_obs         import filter_by_obs_run


app = typer.Typer(no_args_is_help=True)

# --------------------------
# subsetting 
# --------------------------

@app.command("print-dendro-order")
def cmd_print_dendro_order(
    h5ad_in: Path,
    label_key: str = typer.Option(..., "--label-key", "-l"),
):
    """
    Print the dendrogram leaf order (index and label), so you can pick
    positions by number without knowing names.
    """
    order = print_dendro_order_run(h5ad_in, label_key)
    for i, lab in enumerate(order):
        print(f"{i}\t{lab}")

@app.command("subset-by-dendro-range")
def cmd_subset_by_dendro_range(
    h5ad_in: Path,
    h5ad_out: Path,
    label_key: str = typer.Option(..., "--label-key", "-l"),
    start: int = typer.Option(..., "--start", "-s", help="0-based start index (inclusive)"),
    end: int   = typer.Option(..., "--end", "-e", help="0-based end index (inclusive)"),
    invert: bool = typer.Option(False, "--invert", help="Exclude the range instead of keeping it"),
):
    """
    Create a subset .h5ad by selecting a contiguous range of leaves by position
    in the dendrogram (no label strings required).
    """
    subset_by_dendro_range_run(h5ad_in, h5ad_out, label_key, start, end, invert)

# --------------------------
# Plotting
# --------------------------

@app.command("dotplot")
def cmd_dotplot(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    results_csv: Path = typer.Argument(..., exists=True, readable=True),
    label_key: str = typer.Option(..., "--label-key", "-l"),

    # allow column overrides
    cluster_col: str = typer.Option("clusterName", "--cluster-col"),
    markers_col: str = typer.Option("NSForest_markers", "--markers-col", help="Use 'binary_genes' to plot those instead"),

    # outputs
    png_out: Optional[Path] = typer.Option(None, "--png-out"),
    svg_out: Optional[Path] = typer.Option(None, "--svg-out"),
#    html_out: Optional[Path] = typer.Option(None, "--html-out"),
    dpi: int = typer.Option(300, "--dpi"),
):
    if not (png_out or svg_out or html_out):
        raise typer.BadParameter("Provide at least one of --png-out / --svg-out / --html-out.")
    fig = dotplot_run(h5ad_in, results_csv, label_key, cluster_col=cluster_col, markers_col=markers_col)
    if png_out: fig.savefig(str(png_out), bbox_inches="tight", dpi=dpi)
    if svg_out: fig.savefig(str(svg_out), bbox_inches="tight", format="svg")
#    if html_out:
#        pfig = pio.from_matplotlib(fig)
#        pio.write_html(pfig, file=str(html_out), full_html=True, include_plotlyjs="cdn")

@app.command("violinplot")
def cmd_violinplot(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    results_csv: Path = typer.Argument(..., exists=True, readable=True),
    label_key: str = typer.Option(..., "--label-key", "-l"),
    cluster_col: str = typer.Option("clusterName", "--cluster-col"),
    markers_col: str = typer.Option("NSForest_markers", "--markers-col"),

    png_out: Optional[Path] = typer.Option(None, "--png-out"),
    svg_out: Optional[Path] = typer.Option(None, "--svg-out"),
#    html_out: Optional[Path] = typer.Option(None, "--html-out"),
    dpi: int = typer.Option(300, "--dpi"),
):
    if not (png_out or svg_out or html_out):
        raise typer.BadParameter("Provide at least one of --png-out / --svg-out / --html-out.")
    fig = violinplot_run(h5ad_in, results_csv, label_key, cluster_col=cluster_col, markers_col=markers_col)
    if png_out: fig.savefig(str(png_out), bbox_inches="tight", dpi=dpi)
    if svg_out: fig.savefig(str(svg_out), bbox_inches="tight", format="svg")
#    if html_out:
#        pfig = pio.from_matplotlib(fig)
#        pio.write_html(pfig, file=str(html_out), full_html=True, include_plotlyjs="cdn")



@app.command("dendrogramplot")
def cmd_dendrogramplot(
    h5ad_in: Path,
    results_csv: Path,
    h5ad_out: Path,
    label_key: str = typer.Option(..., "--label-key", "-l"),
    png_out: Optional[Path] = typer.Option(None, "--png-out", help="Save dendrogram as PNG"),
    svg_out: Optional[Path] = typer.Option(None, "--svg-out", help="Save dendrogram as SVG"),
):
    adata = dendrogramplot_run(
        h5ad_in, results_csv, label_key,
        png_out=png_out, svg_out=svg_out
    )
    if adata is not None:
        adata.write_h5ad(str(h5ad_out))

# --------------------------
# Preprocessing helpers
# --------------------------

@app.command("prep-medians")
def cmd_prep_medians(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    h5ad_out: Path = typer.Argument(...),
    label_key: str = typer.Option(..., "--label-key", "-l"),
):
    """
    Compute per-cluster medians (nsforest.pp.prep_medians) and write a new .h5ad.
    """
    adata = prep_medians_run(h5ad_in, label_key)
    adata.write_h5ad(str(h5ad_out))


@app.command("prep-binary-scores")
def cmd_prep_binary_scores(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    h5ad_out: Path = typer.Argument(...),
    label_key: str = typer.Option(..., "--label-key", "-l"),
):
    """
    Compute per-cluster binary scores (nsforest.pp.prep_binary_scores) and write a new .h5ad.
    """
    adata = prep_binary_scores_run(h5ad_in, label_key)
    adata.write_h5ad(str(h5ad_out))


# --------------------------
# Core NS-Forest
# --------------------------

@app.command("nsforest")
def cmd_nsforest(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    results_csv: Path = typer.Argument(...),
    label_key: str = typer.Option(..., "--label-key", "-l"),
    # passthrough tuning params (mirrors upstream)
    output_folder: Optional[Path] = typer.Option(None, "--output-folder", help="If provided, NS-Forest may write internals there."),
    cluster_list: Optional[str] = typer.Option(None, "--cluster-list", help='JSON list or comma-separated.'),
    n_trees: int = typer.Option(1000, "--n-trees"),
    n_jobs: int = typer.Option(-1, "--n-jobs"),
    beta: float = typer.Option(0.5, "--beta"),
    n_top_genes: int = typer.Option(15, "--n-top-genes"),
    n_binary_genes: int = typer.Option(10, "--n-binary-genes"),
    n_genes_eval: int = typer.Option(6, "--n-genes-eval"),
    save_supplementary: bool = typer.Option(False, "--save-supplementary"),
):
    """
    Run NS-Forest core algorithm and write the results CSV.
    This CLI does not create directories; if output_folder is used, ensure it exists upstream.
    """
    df = nsforest_run(
        h5ad_in=h5ad_in,
        label_key=label_key,
        output_folder=output_folder,
        cluster_list=cluster_list,
        n_trees=n_trees,
        n_jobs=n_jobs,
        beta=beta,
        n_top_genes=n_top_genes,
        n_binary_genes=n_binary_genes,
        n_genes_eval=n_genes_eval,
        save_supplementary=save_supplementary,
    )
    df.to_csv(results_csv, index=False)


# --------------------------
# Evaluation
# --------------------------

@app.command("eval-markers")
def cmd_eval_markers(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    markers_json: Path = typer.Argument(..., exists=True, readable=True),
    results_csv: Path = typer.Argument(...),
    label_key: str = typer.Option(..., "--label-key", "-l"),
):
    """
    Evaluate marker sets using nsforest.evaluating helpers; write results CSV.
    """
    df = eval_markers_run(h5ad_in, markers_json, label_key)
    df.to_csv(results_csv, index=False)


# --------------------------
# Utilities: sanitize & filter
# --------------------------

@app.command("sanitize-labels")
def cmd_sanitize_labels(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    h5ad_out: Path = typer.Argument(...),
    label_key: str = typer.Option(..., "--label-key", "-l"),
):
    """
    Collision-safe label sanitizer (replaces non [A-Za-z0-9_-] with '_', collapses repeats, dedupes).
    Writes a new .h5ad to the path you provide.
    """
    adata = sanitize_labels_run(h5ad_in, label_key)
    adata.write_h5ad(str(h5ad_out))


@app.command("filter-by-obs")
def cmd_filter_by_obs(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    h5ad_out: Path = typer.Argument(...),
    obs_key: str = typer.Option(..., "--obs-key", "-k", help="obs column to filter on (e.g., 'tissue' or 'disease')"),
    values: List[str] = typer.Option(..., "--value", "-v", help="Value(s) to KEEP (repeat flag for multiples)"),
    mode: str = typer.Option("exact", "--mode", help="exact | contains | regex"),
    case_insensitive: bool = typer.Option(True, "--ci/--no-ci"),
    na_policy: str = typer.Option("drop", "--na", help="drop | keep | match"),
    invert: bool = typer.Option(False, "--invert", help="Keep non-matching rows if True"),
):
    """
    Filter by a SINGLE obs field. Chain this command multiple times for tissue, then disease, etc.
    Writes the filtered .h5ad to the path you provide.
    """
    adata = filter_by_obs_run(
        h5ad_in=h5ad_in,
        obs_key=obs_key,
        values=values,
        mode=mode, case_insensitive=case_insensitive,
        na_policy=na_policy, invert=invert,
    )
    adata.write_h5ad(str(h5ad_out))


# --------------------------
# Entrypoint
# --------------------------

if __name__ == "__main__":
    app()

