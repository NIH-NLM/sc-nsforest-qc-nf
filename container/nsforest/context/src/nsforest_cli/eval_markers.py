from pathlib import Path
import json
import pandas as pd
import scanpy as sc
from nsforest import ns, nsforesting, utils, NSFOREST_VERSION


def eval_markers_run(
        *,
        h5ad_in: Path,
        markers_csv: Path,
        label_key: str,

) -> pd.DataFrame:
    """
    Evaluate user marker sets using nsforest.evaluating helpers.
    Returns a DataFrame; caller writes CSV.
    """
    adata = sc.read_h5ad(str(h5ad_in))
    markers = pd.read_csv(str(markers_csv))
    markers_dict = utils.prepare_markers(markers, "clusterName", "markers")
    cluster_header        = "clusterName"
    outputfilename_prefix = "marker_eval"
    output_folder         = "."
    evaluation_results = ns.ev.DecisionTree(
        adata,
        cluster_header,
        markers_dict,
        combinations = True,
        use_mean = False,
        save = True,
        save_supplementary = True,
        output_folder = output_folder,
        outputfilename_prefix = outputfilename_prefix,
    )

    return evaluation_results

