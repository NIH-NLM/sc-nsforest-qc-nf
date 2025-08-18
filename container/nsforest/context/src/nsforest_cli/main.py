# container/nsforest/context/src/nsforest_cli/main.py
from __future__ import annotations

from pathlib import Path
from typing import List, Optional
import typer
import plotly.io as pio

# ---- pure helpers: each provides exactly one function ----
from nsforest_cli.dotplot import dotplot_run
from nsforest_cli.violinplot import violinplot_run
from nsforest_cli.dendrogramplot import dendrogramplot_run

from nsforest_cli.prep_medians import prep_medians_run
from nsforest_cli.prep_binary_scores import prep_binary_scores_run
from nsforest_cli.run_nsforest import nsforest_run
from nsforest_cli.eval_markers import eval_markers_run

from nsforest_cli.sanitize_labels import sanitize_labels_run
from nsforest_cli.filter_obs import filter_by_obs_run

app = typer.Typer(no_args_is_help=True)


# --------------------------
# Plotting
# --------------------------

@app.command("dotplot")
def cmd_dotplot(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    markers_json: Path = typer.Argument(..., exists=True, readable=True),
    label_key: str = typer.Option(..., "--label-key", "-l", help="obs column with cluster labels"),
    # outputs: require caller to specify at least one
    png_out: Optional[Path] = typer.Option(None, "--png-out", help="Write static PNG via Matplotlib"),
    svg_out: Optional[Path] = typer.Option(None, "--svg-out", help="Write static SVG via Matplotlib"),
    html_out: Optional[Path] = typer.Option(None, "--html-out", help="Write interactive HTML via Plotly"),
    dpi: int = typer.Option(300, "--dpi"),
):
    """
    Dotplot using NS-Forest plotting helper.
    Writes only the files you explicitly request. No directories are created.
    """
    if not (png_out or svg_out or html_out):
        raise typer.BadParameter("Provide at least one of --png-out / --svg-out / --html-out.")

    fig = dotplot_run(h5ad_in, markers_json, label_key)

    if png_out:
        fig.savefig(str(png_out), bbox_inches="tight", dpi=dpi)
    if svg_out:
        fig.savefig(str(svg_out), bbox_inches="tight", format="svg")
    if html_out:
        pfig = pio.from_matplotlib(fig)
        pio.write_html(pfig, file=str(html_out), full_html=True, include_plotlyjs="cdn")


@app.command("violinplot")
def cmd_violinplot(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    markers_json: Path = typer.Argument(..., exists=True, readable=True),
    label_key: str = typer.Option(..., "--label-key", "-l"),
    png_out: Optional[Path] = typer.Option(None, "--png-out"),
    svg_out: Optional[Path] = typer.Option(None, "--svg-out"),
    html_out: Optional[Path] = typer.Option(None, "--html-out"),
    dpi: int = typer.Option(300, "--dpi"),
):
    """
    Stacked violin using NS-Forest plotting helper.
    Writes only the files you explicitly request. No directories are created.
    """
    if not (png_out or svg_out or html_out):
        raise typer.BadParameter("Provide at least one of --png-out / --svg-out / --html-out.")

    fig = violinplot_run(h5ad_in, markers_json, label_key)

    if png_out:
        fig.savefig(str(png_out), bbox_inches="tight", dpi=dpi)
    if svg_out:
        fig.savefig(str(svg_out), bbox_inches="tight", format="svg")
    if html_out:
        pfig = pio.from_matplotlib(fig)
        pio.write_html(pfig, file=str(html_out), full_html=True, include_plotlyjs="cdn")


# --------------------------
# Dendrogram metadata
# --------------------------

@app.command("dendrogramplot")
def cmd_dendrogramplot(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    h5ad_out: Path = typer.Argument(...),
    label_key: str = typer.Option(..., "--label-key", "-l"),
):
    """
    Compute/attach dendrogram metadata via NS-Forest preprocessing helper.
    Writes the updated .h5ad to the path you provide.
    """
    adata = dendrogramplot_run(h5ad_in, label_key)
    adata.write_h5ad(str(h5ad_out))


# --------------------------
# Preprocessing helpers
# --------------------------

@app.command("prep-medians")
def cmd_prep_medians(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    h5ad_out: Path = typer.Argument(...),
    label_key: str = typer.Option(..., "--label-key", "-l"),
    medians_header: str = typer.Option("medians_", "--medians-header"),
):
    """
    Compute per-cluster medians (nsforest.pp.prep_medians) and write a new .h5ad.
    """
    adata = prep_medians_run(h5ad_in, label_key, medians_header=medians_header)
    adata.write_h5ad(str(h5ad_out))


@app.command("prep-binary-scores")
def cmd_prep_binary_scores(
    h5ad_in: Path = typer.Argument(..., exists=True, readable=True),
    h5ad_out: Path = typer.Argument(...),
    label_key: str = typer.Option(..., "--label-key", "-l"),
    binary_scores_header: str = typer.Option("binary_scores_", "--binary-scores-header"),
):
    """
    Compute per-cluster binary scores (nsforest.pp.prep_binary_scores) and write a new .h5ad.
    """
    adata = prep_binary_scores_run(h5ad_in, label_key, binary_scores_header=binary_scores_header)
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

