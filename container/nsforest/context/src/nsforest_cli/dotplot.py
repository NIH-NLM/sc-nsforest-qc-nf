from pathlib import Path
import scanpy as sc
import nsforest as ns
from nsforest_cli._markers_from_results import markers_from_nsforest_results

def _as_figure(obj):
    # Try to extract a Matplotlib Figure from whatever was returned
    if obj is None:
        return None
    # scanpy’s DotPlot object
    if hasattr(obj, "fig"):
        return obj.fig
    # some wrappers hand back a matplotlib Axes
    if hasattr(obj, "figure"):
        return obj.figure
    # some scanpy plotters return a matplotlib Figure directly
    try:
        import matplotlib.figure
        if isinstance(obj, matplotlib.figure.Figure):
            return obj
    except Exception:
        pass
    return None

def dotplot_run(
    h5ad_in: Path,
    results_csv: Path,
    label_key: str,
    cluster_col: str = "clusterName",
    markers_col: str = "NSForest_markers",  # or "binary_genes"
):
    """
    Return a Matplotlib Figure for dotplot using NS-Forest plotting,
    deriving markers directly from results.csv (no JSON required).
    Robust to differing return types from scanpy/nsforest.
    """
    adata = sc.read_h5ad(str(h5ad_in))
    markers = markers_from_nsforest_results(results_csv, cluster_col=cluster_col, markers_col=markers_col)

    # First try NS-Forest’s wrapper
    obj = ns.plotting.dotplot(adata, markers, label_key)
    fig = _as_figure(obj)
    if fig is not None:
        return fig

    # Fallback: call scanpy directly to ensure we get a handle
    dp = sc.pl.dotplot(adata, markers, groupby=label_key, show=False)  # returns a DotPlot
    # DotPlot has .fig in modern scanpy; if not, grab current figure
    fig = getattr(dp, "fig", None)
    if fig is None:
        import matplotlib.pyplot as plt
        fig = plt.gcf()
    return fig

