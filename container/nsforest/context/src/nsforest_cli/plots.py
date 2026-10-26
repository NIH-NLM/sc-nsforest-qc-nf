"""
Create NSForest visualization plots.

Corresponds to DEMO_NS-Forest_workflow.py: Section 4

Saves boxplots (html), scatter plots (svg), and expression plots (svg).
"""

import matplotlib
matplotlib.use("Agg")

import ast
import glob
import os
import scanpy as sc
import pandas as pd
from . import nsf_compat

from .common_utils import (
    get_output_prefix,
    load_h5ad,
    log_section,
    logger
)
from .gene_mapping_utils import (
    load_gene_mapping,
    map_markers_to_symbols,
    add_gene_symbols_to_adata
)


def run_plots(h5ad_path, results_csv, cluster_header, organ, first_author, journal, year, embedding, dataset_version_id,
              max_cells_per_cluster=0, seed=42):
    """
    Create NSForest visualization plots with gene symbol mapping.
    """
    log_section("NSForest: Plotting")
    sc.settings.figdir = "."

    prefix = get_output_prefix( organ, first_author, journal, year, cluster_header, embedding, dataset_version_id )

    # Load results
    logger.info(f"Loading results: {results_csv}")

    results = pd.read_csv(results_csv)
    results = results.dropna(subset=['NSForest_markers', 'clusterName'])
    results = results[results['NSForest_markers'].str.strip() != '[]']

    if results.empty:
        logger.warning("No valid marker results — skipping all plots")
        return

    results['NSForest_markers'] = results['NSForest_markers'].apply(ast.literal_eval)
    
    # Gene symbol mapping
    ensg_to_symbol = load_gene_mapping()
    results, markers_dict = map_markers_to_symbols(results, ensg_to_symbol)

    # Boxplots — html (interactive) + svg (static for publication)
    for metric in ['f_score', 'precision', 'recall', 'onTarget']:
        nsf_compat.boxplot(results, metric, save=True, output_folder="", outputfilename_prefix=prefix)
        os.rename(f"{prefix}_boxplot_{metric}.html", f"boxplot_{prefix}_{metric}.html")
        nsf_compat.boxplot(results, metric, save="svg",  output_folder="", outputfilename_prefix=prefix)
        os.rename(f"{prefix}_boxplot_{metric}.svg",  f"boxplot_{prefix}_{metric}.svg")
    logger.info("Boxplots saved.")

    # Scatter plots — html (interactive) + svg (static for publication)
    for metric in ['f_score', 'precision', 'recall', 'onTarget']:
        nsf_compat.scatter_w_clusterSize(results, metric, save=True,  output_folder="", outputfilename_prefix=prefix)
        os.rename(f"{prefix}_scatter_{metric}.html", f"scatter_{prefix}_{metric}.html")
        nsf_compat.scatter_w_clusterSize(results, metric, save="svg", output_folder="", outputfilename_prefix=prefix)
        os.rename(f"{prefix}_scatter_{metric}.svg",  f"scatter_{prefix}_{metric}.svg")
    logger.info("Scatter plots saved.")
    
    # Load adata for expression plots
    logger.info(f"Loading h5ad: {h5ad_path}")
    adata = load_h5ad(h5ad_path, cluster_header)
    adata = add_gene_symbols_to_adata(adata, ensg_to_symbol)

    # some cells in adata.obs[cluster_header] have NaN (float) instead of a string label.
    # The dendrogram reorder will fail because it cant join floats as strings...
    adata = adata[adata.obs[cluster_header].notna()].copy()
    adata.obs[cluster_header] = adata.obs[cluster_header].astype(str).astype("category")

    # Cap cells per cluster for the expression plots — dotplot/violin/matrix aggregate per cluster,
    # so a few thousand cells/cluster is plenty, while the full multi-million-cell matrix OOMs
    # plotting. Stratified + seeded; clusters smaller than the cap keep all their cells.
    if max_cells_per_cluster and max_cells_per_cluster > 0:
        import numpy as np
        rng = np.random.default_rng(seed)
        groups = adata.obs.groupby(cluster_header, observed=True).indices
        keep = []
        for cl, idx in groups.items():
            idx = np.asarray(idx)
            if len(idx) > max_cells_per_cluster:
                idx = rng.choice(idx, size=max_cells_per_cluster, replace=False)
            keep.append(idx)
        keep_idx = np.sort(np.concatenate(keep))
        n_before = adata.n_obs
        adata = adata[keep_idx].copy()
        logger.info(f"Subsampled cells for plotting: {n_before} -> {adata.n_obs} "
                    f"(max {max_cells_per_cluster}/cluster, seed={seed})")

    # Dotplot
    nsf_compat.dotplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                  gene_symbols='gene_symbol', save="svg", output_folder="",
                  outputfilename_suffix=prefix)
    nsf_compat.dotplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                  gene_symbols='gene_symbol', standard_scale='var', save="svg",
                  output_folder="", outputfilename_suffix="_scaled" + prefix)

    # Stacked violin
    nsf_compat.stackedviolin(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                        gene_symbols='gene_symbol', save="svg", output_folder="",
                        outputfilename_suffix=prefix)
    nsf_compat.stackedviolin(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                        gene_symbols='gene_symbol', standard_scale='var', save="svg",
                        output_folder="", outputfilename_suffix="_scaled" + prefix)

    # Matrix plot
    nsf_compat.matrixplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                     gene_symbols='gene_symbol', save="svg", output_folder="",
                     outputfilename_suffix=prefix)
    nsf_compat.matrixplot(adata, markers_dict, cluster_header, dendrogram=True, use_raw=False,
                     gene_symbols='gene_symbol', standard_scale='var', save="svg",
                     output_folder="", outputfilename_suffix="_scaled" + prefix)


    logger.info("Plotting complete!")


def parent_group_order(obs, parent_header, child_header, dend_order_csv, order_column="CrossArea_cluster_order"):
    """
    The order of the parent groups (the subclasses) of a figure of the paper: the order in which the clusters
    of the dendrogram file (Jorstad Fig 5B) first meet each parent.
    """
    parents = obs[[child_header, parent_header]].astype(str).drop_duplicates()
    # the dendrogram file may spell a cluster in capitals (PAX6_1 for Pax6_1): the names are matched without case
    parent_of = {child.casefold(): parent for child, parent in zip(parents[child_header], parents[parent_header])}
    order = []
    for child in pd.read_csv(dend_order_csv)[order_column].astype(str):
        parent = parent_of.get(child.casefold())
        if parent is not None and parent not in order:
            order.append(parent)
    return order


def _order_from_file(names, dend_order_csv, order_column="CrossArea_cluster_order"):
    """`names` in the order of the dendrogram file (matched without case: PAX6_1 is Pax6_1); the names it does not mention follow."""
    by_fold = {n.casefold(): n for n in names}
    ordered = []
    for n in pd.read_csv(dend_order_csv)[order_column].astype(str):
        real = by_fold.get(n.casefold())
        if real is not None and real not in ordered:
            ordered.append(real)
    return ordered + [n for n in names if n not in ordered]


def _load_for_global_plots(h5ad_path, results_csv, cluster_header):
    """The markers of the results (as gene symbols) and the h5ad with its gene symbols, for the plots of the paper."""
    logger.info(f"Loading results: {results_csv}")
    results = pd.read_csv(results_csv)
    results = results.dropna(subset=['NSForest_markers', 'clusterName'])
    results = results[results['NSForest_markers'].str.strip() != '[]']
    if results.empty:
        return None, None
    results['NSForest_markers'] = results['NSForest_markers'].apply(ast.literal_eval)

    ensg_to_symbol = load_gene_mapping()
    results, markers_dict = map_markers_to_symbols(results, ensg_to_symbol)

    logger.info(f"Loading h5ad: {h5ad_path}")
    adata = load_h5ad(h5ad_path, cluster_header)
    adata = add_gene_symbols_to_adata(adata, ensg_to_symbol)
    adata = adata[adata.obs[cluster_header].notna()].copy()
    adata.obs[cluster_header] = adata.obs[cluster_header].astype(str).astype("category")
    return adata, markers_dict


def _marker_labels(plot):
    """For each column of a scanpy plot: the gene and the group of markers (the row it belongs to) it stands for."""
    genes = [str(g) for g in plot.dot_size_df.columns] if hasattr(plot, "dot_size_df") else [str(g) for g in plot.values_df.columns]
    groups = [""] * len(genes)
    labels, positions = getattr(plot, "var_group_labels", None), getattr(plot, "var_group_positions", None)
    if labels is not None and positions is not None:
        for label, (start, end) in zip(labels, positions):
            for k in range(start, min(end, len(genes) - 1) + 1):
                groups[k] = str(label)
    return genes, groups


def _html_dotplot(plot, path, title, row_order=None):
    """Interactive dotplot: size = fraction of cells expressing the marker, colour = scaled mean expression."""
    import numpy as np
    import plotly.graph_objects as go
    size, color = plot.dot_size_df, plot.dot_color_df
    if row_order is not None:   # scanpy applies categories_order only when it draws
        size, color = size.loc[row_order], color.loc[row_order]
    rows = [str(r) for r in size.index]
    genes, marker_of = _marker_labels(plot)
    n_r, n_c = size.shape
    frac, val = size.values.ravel(), color.values.ravel()
    custom = np.column_stack([np.repeat(rows, n_c), np.tile(genes, n_r), np.tile(marker_of, n_r),
                              np.round(frac, 3).astype(str), np.round(val, 3).astype(str)])
    fig = go.Figure(go.Scattergl(
        x=np.tile(np.arange(n_c), n_r), y=np.repeat(np.arange(n_r), n_c), mode="markers",
        marker=dict(size=3 + 12 * frac, color=val, colorscale="Reds", cmin=0, cmax=1,
                    colorbar=dict(title="scaled mean<br>expression"), line=dict(width=0)),
        customdata=custom,
        hovertemplate=("%{customdata[0]}<br>marker %{customdata[1]} (of %{customdata[2]})<br>"
                       "fraction of cells %{customdata[3]}<br>scaled mean %{customdata[4]}<extra></extra>")))
    fig.update_layout(title=title, width=max(900, 7 * n_c), height=max(500, 14 * n_r), plot_bgcolor="white")
    fig.update_yaxes(tickmode="array", tickvals=list(range(n_r)), ticktext=rows, autorange="reversed")
    fig.update_xaxes(tickmode="array", tickvals=list(range(n_c)), ticktext=genes, tickangle=90, side="top")
    fig.write_html(path, include_plotlyjs="cdn")


def _html_matrixplot(plot, path, title, row_order=None):
    """Interactive matrix plot: scaled mean expression of each marker in each row."""
    import numpy as np
    import plotly.graph_objects as go
    values = plot.values_df
    if row_order is not None:
        values = values.loc[row_order]
    rows = [str(r) for r in values.index]
    genes, marker_of = _marker_labels(plot)
    n_r, n_c = values.shape
    text = np.array([[f"{g} (of {m})" for g, m in zip(genes, marker_of)]] * n_r)
    fig = go.Figure(go.Heatmap(z=np.round(values.values, 3), x=list(range(n_c)), y=rows, text=text, colorscale="Viridis",
                               zmin=0, zmax=1, colorbar=dict(title="scaled mean<br>expression"),
                               hovertemplate="%{y}<br>marker %{text}<br>scaled mean %{z}<extra></extra>"))
    fig.update_layout(title=title, width=max(900, 7 * n_c), height=max(500, 14 * n_r))
    fig.update_yaxes(autorange="reversed")
    fig.update_xaxes(tickmode="array", tickvals=list(range(n_c)), ticktext=genes, tickangle=90, side="top")
    fig.write_html(path, include_plotlyjs="cdn")


def _html_picture(base, title, use_svg):
    """An HTML page that holds the picture of a figure: the SVG inline, or the PNG when the SVG is too large to embed."""
    import base64
    if use_svg:
        svg = open(f"{base}.svg", encoding="utf-8").read()
        svg = svg[svg.index("<svg"):]
        body = svg
    else:
        data = base64.b64encode(open(f"{base}.png", "rb").read()).decode()
        body = f'<img src="data:image/png;base64,{data}" alt="{title}">'
    with open(f"{base}.html", "w", encoding="utf-8") as fh:
        fh.write(f"<!doctype html><html><head><meta charset='utf-8'><title>{title}</title></head>"
                 f"<body><h3>{title}</h3>{body}</body></html>")


def _draw_global_figures(adata, markers_dict, groupby, order, prefix, stem, svg_in_html):
    """
    The dotplot, the stacked violin and the matrix plot, each with every marker scaled to 0-1 across the rows
    (standard_scale='var'; "scaled" is in the file names), each as PNG (100 dpi, as the authors' PNG files), SVG and HTML:
    <plot>_<stem>_scaled_<prefix>.{png,svg,html}. A figure is drawn once and saved in both picture formats. The HTML
    of the dotplot and of the matrix plot is interactive (plotly), the HTML of the stacked violin holds the picture.
    """
    import matplotlib.pyplot as plt
    kinds = (("dotplot", sc.pl.dotplot), ("stacked_violin", sc.pl.stacked_violin), ("matrixplot", sc.pl.matrixplot))
    row_order = dict(dendrogram=True) if order is True else dict(dendrogram=False, categories_order=order)
    for kind, draw in kinds:
        base = f"{kind}_{stem}_scaled_{prefix}"
        title = f"{kind} {stem} (scaled): {prefix}"
        logger.info(f"Drawing {base}")
        plot = draw(adata, markers_dict, groupby, use_raw=False, gene_symbols='gene_symbol', standard_scale='var',
                    return_fig=True, show=False, **row_order)
        if kind == "dotplot":
            _html_dotplot(plot, f"{base}.html", title, None if order is True else order)
        elif kind == "matrixplot":
            _html_matrixplot(plot, f"{base}.html", title, None if order is True else order)
        plot.make_figure()
        plt.savefig(f"{base}.png", bbox_inches="tight", dpi=100)
        plt.savefig(f"{base}.svg", bbox_inches="tight")
        plt.close("all")
        if kind == "stacked_violin":
            _html_picture(base, title, use_svg=svg_in_html)


def run_subclass_plots(h5ad_path, results_csv, cluster_header, child_header, organ, first_author, journal, year, embedding,
                       dataset_version_id, dend_order_csv=None, order_column="CrossArea_cluster_order"):
    """
    The global subclass figures of the Zenodo deposit (stacked_violin_global_subclass.png, dotplot_global_subclass.png)
    and the matrix plot: one row per subclass, and under each subclass its own NS-Forest markers, on all the cells.
    The rows follow the order of the dendrogram file when it is given.
    """
    log_section("NSForest: Plotting the global subclass results (scaled dotplot, stacked violin, matrix plot)")
    prefix = get_output_prefix(organ, first_author, journal, year, cluster_header, embedding, dataset_version_id)
    adata, markers_dict = _load_for_global_plots(h5ad_path, results_csv, cluster_header)
    if adata is None:
        logger.warning("No valid marker results: skipping the subclass plots")
        return

    order = True
    if dend_order_csv:
        order = parent_group_order(adata.obs, cluster_header, child_header, dend_order_csv, order_column)
        order = [g for g in order if g in markers_dict] + [g for g in markers_dict if g not in order]
        markers_dict = {g: markers_dict[g] for g in order}
        logger.info(f"Subclass rows in the order of {dend_order_csv}: {order}")

    _draw_global_figures(adata, markers_dict, cluster_header, order, prefix, "global_subclass", svg_in_html=True)
    logger.info("Subclass plotting complete!")


def run_combinatorial_plots(h5ad_path, results_csv, cluster_header, organ, first_author, journal, year, embedding,
                            dataset_version_id, dend_order_csv=None, order_column="CrossArea_cluster_order"):
    """
    The global cluster combinatorial figures (dotplot_global_cluster_combinatorial_on_global.png of the Zenodo deposit,
    and the stacked violin and the matrix plot): one row per cluster, and under each cluster its best combination of
    subclass and cluster markers, on all the cells. The rows follow the order of the dendrogram file when it is given.
    The SVG of the stacked violin of 153 clusters is hundreds of MB: its HTML page holds the PNG.
    """
    log_section("NSForest: Plotting the global combinatorial results (scaled dotplot, stacked violin, matrix plot)")
    prefix = get_output_prefix(organ, first_author, journal, year, cluster_header, embedding, dataset_version_id)
    adata, markers_dict = _load_for_global_plots(h5ad_path, results_csv, cluster_header)
    if adata is None:
        logger.warning("No valid marker results: skipping the combinatorial plots")
        return

    order = True
    if dend_order_csv:
        order = _order_from_file(list(markers_dict), dend_order_csv, order_column)
        markers_dict = {g: markers_dict[g] for g in order}
        logger.info(f"Cluster rows in the order of {dend_order_csv}: {len(order)} clusters")

    _draw_global_figures(adata, markers_dict, cluster_header, order, prefix, "global_cluster_combinatorial_on_global",
                         svg_in_html=False)
    logger.info("Combinatorial plotting complete!")
