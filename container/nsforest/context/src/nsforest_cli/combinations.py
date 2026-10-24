"""
Evaluate marker combinations per cluster (the combinatorial step of Jorstad et al. 2023).

Corresponds to the Zenodo notebook "NS-Forest markers for human neocortex cross-area cluster":
the NS-Forest markers of the subclasses and of the clusters are combined and evaluated with
nsforest.evaluating.DecisionTree(..., combinations=True), which for each cluster keeps the
combination of its candidate genes with the best F-beta score.

A cluster's candidate genes are its own NS-Forest markers followed by the NS-Forest markers of the
subclass it belongs to (the parent is read from adata.obs). Genes appear once, in that order.

Saves:
  results_combinations_{organ}_{first_author}_{journal}_{year}_{cluster_header}_{embedding}_{vid}.csv
"""
import ast
import inspect
import pandas as pd

from .common_utils import get_output_prefix, load_h5ad, log_section, logger


def _read_markers(results_csv):
    """clusterName -> list of NSForest markers, from a merged NSForest results CSV."""
    df = pd.read_csv(results_csv)
    markers = {}
    for name, value in zip(df["clusterName"], df["NSForest_markers"]):
        markers[str(name)] = ast.literal_eval(value) if isinstance(value, str) else list(value)
    return markers


def build_markers_dict(adata, cluster_header, subclass_header, cluster_markers, subclass_markers):
    """Candidate genes per cluster: its own markers, then its parent subclass's, each gene once."""
    parents = (
        adata.obs[[cluster_header, subclass_header]].astype(str).drop_duplicates()
    )
    multi = parents.groupby(cluster_header)[subclass_header].nunique()
    if (multi > 1).any():
        raise ValueError(f"Clusters with more than one parent {subclass_header}: {list(multi[multi > 1].index)}")
    parent_of = dict(zip(parents[cluster_header], parents[subclass_header]))

    markers_dict = {}
    for cluster in sorted(adata.obs[cluster_header].astype(str).unique()):
        genes = list(cluster_markers.get(cluster, []))
        for gene in subclass_markers.get(parent_of[cluster], []):
            if gene not in genes:
                genes.append(gene)
        if genes:
            markers_dict[cluster] = genes
        else:
            logger.warning(f"No candidate markers for cluster {cluster}; it is left out")
    return markers_dict


def run_combinations(h5ad_path, cluster_header, subclass_header, subclass_results_csv, cluster_results_csv,
                     organ, first_author, journal, year, embedding, dataset_version_id):
    """Run DecisionTree with combinations=True over the combined subclass and cluster markers."""
    log_section("NSForest: marker combinations (DecisionTree, combinations=True)")
    from importlib.metadata import version as _pkg_version
    from nsforest import evaluating as ev
    logger.info(f"NS-Forest version: {_pkg_version('nsforest')}")

    prefix = get_output_prefix(organ, first_author, journal, year, cluster_header, embedding, dataset_version_id)

    adata = load_h5ad(h5ad_path, cluster_header)
    if subclass_header not in adata.obs.columns:
        raise ValueError(f"Subclass column '{subclass_header}' not found in adata.obs")

    subclass_markers = _read_markers(subclass_results_csv)
    cluster_markers = _read_markers(cluster_results_csv)
    logger.info(f"Markers read: {len(subclass_markers)} subclasses, {len(cluster_markers)} clusters")

    markers_dict = build_markers_dict(adata, cluster_header, subclass_header, cluster_markers, subclass_markers)
    logger.info(f"Evaluating combinations for {len(markers_dict)} clusters")

    kwargs = dict(combinations=True)
    if 'save' in inspect.signature(ev.DecisionTree).parameters:
        kwargs['save'] = False
    else:
        # 4.0.0 always creates output_folder and writes into it; keep that out of the task outputs
        kwargs['output_folder'] = "nsforest_native_output/"
        kwargs['outputfilename_prefix'] = "native"

    results = ev.DecisionTree(adata, cluster_header, markers_dict, **kwargs)
    logger.info(f"Combinations results shape: {results.shape}")

    output_csv = f"results_combinations_{prefix}.csv"
    results.to_csv(output_csv, index=False)
    logger.info(f"Saved: {output_csv}")
    logger.info("Combinations complete!")
