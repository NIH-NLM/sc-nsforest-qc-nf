from pathlib import Path
from typing import List, Optional
import json
import pandas as pd
import scanpy as sc
import nsforest as ns

def nsforest_run(
        *,
        h5ad_in: Path,
        label_key: str,
        results_csv_out: Path,
        output_folder: Optional[Path] = None,
        cluster_list: Optional[str] = None,
        n_trees: int = 1000,
        n_jobs: int = -1,
        beta: float = 0.5,
        n_top_genes: int = 15,
        n_binary_genes: int = 10,
        n_genes_eval: int = 6,
        save_supplementary: bool = True,
) -> pd.DataFrame:
    """
    Run NS-Forest core algorithm and return the results DataFrame.
    No I/O here; caller decides where to write CSV. Pass output_folder if upstream expects it.
    """
    adata = sc.read_h5ad(str(h5ad_in))
    cl: List[str] = []

    if cluster_list:
        cl = json.loads(cluster_list) if cluster_list.strip().startswith("[") else [x for x in cluster_list.split(",") if x]

    # Set a default output folder if not provided
    if output_folder is None:
        output_folder = Path("nsforest_outputs")
    output_folder.mkdir(parents=True, exist_ok=True)

    df = ns.nsforesting.NSForest(
        adata,
        cluster_header=label_key,
        cluster_list=cl,
        n_trees=n_trees,
        n_jobs=n_jobs,
        beta=beta,
        n_top_genes=n_top_genes,
        n_binary_genes=n_binary_genes,
        n_genes_eval=n_genes_eval,
        save=True,
        save_supplementary=save_supplementary,
        output_folder=str(output_folder),
    )

    df.to_csv(results_csv_out, index=False)
    return None

