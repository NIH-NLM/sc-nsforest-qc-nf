from pathlib import Path
import scanpy as sc
import nsforest as ns
from nsforest_cli._markers_from_results import markers_from_nsforest_results

def _as_figure(obj):
    if obj is None:
        return None
    if hasattr(obj, "fig"):
        return obj.fig
    if hasattr(obj, "figure"):
        return obj.figure
    try:
        import matplotlib.figure
        if isinstance(obj, matplotlib.figure.Figure):
            return obj
    except Exception:
        pass
    return None

def violinplot_run(
    h5ad_in: Path,
    results_csv: Path,
    label_key: str,
    cluster_col: str = "clusterName",
    markers_col: str = "NSForest_markers",  # or "binary_genes"
):
    """
    Return a Matplotlib Figure for stacked violin using NS-Forest plotting,
    deriving markers directly from results.csv. Robust to return-type differences.
    """
    adata = sc.read_h5ad(str(h5ad_in))
    markers = markers_from_nsforest_results(results_csv, cluster_col=cluster_col, markers_col=markers_col)

    obj = ns.plotting.stackedviolin(adata, markers, label_key)
    fig = _as_figure(obj)
    if fig is not None:
        return fig

    # Fallback to scanpy
    vg = sc.pl.stacked_violin(adata, markers, groupby=label_key, show=False)  # returns a matplotlib Figure in newer versions
    # Newer scanpy returns a Figure; if older, take current figure
    try:
        import matplotlib.figure
        if isinstance(vg, matplotlib.figure.Figure):
            return vg
    except Exception:
        pass

    import matplotlib.pyplot as plt
    return plt.gcf()

