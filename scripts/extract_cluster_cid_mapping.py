"""Extract a per-cluster cell-ontology mapping table from an .h5ad file.

Reads an AnnData file and, for each unique value in the cluster-header column
of ``adata.obs``, emits a row pairing that cluster with the cell-ontology term
currently assigned in the AnnData. Output is a CSV with exactly four columns:

    cluster_name, skos, manual_mapped_cid, cell_ontology_id

- ``cluster_name``       : unique value from the cluster-header column.
- ``skos``               : "exact" or "related" (restricted vocabulary).
- ``manual_mapped_cid``  : blank, to be filled in manually downstream.
- ``cell_ontology_id``   : value from the cell-ontology column in obs
                           (defaults to ``cell_type_ontology_term_id``).

If a cluster is associated with more than one cell_ontology_id in obs, one
row is emitted per (cluster, cid) pair.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc

ALLOWED_SKOS = {"exact", "related"}


def build_mapping(
    h5ad_path: str,
    cluster_header: str,
    cid_column: str = "cell_type_ontology_term_id",
    default_skos: str = "",
) -> pd.DataFrame:
    if default_skos and default_skos not in ALLOWED_SKOS:
        raise ValueError(
            f"--default-skos must be one of {sorted(ALLOWED_SKOS)}, got {default_skos!r}"
        )

    adata = sc.read_h5ad(h5ad_path)

    if cluster_header not in adata.obs.columns:
        raise KeyError(
            f"cluster-header '{cluster_header}' not found in obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )
    if cid_column not in adata.obs.columns:
        raise KeyError(
            f"cell-ontology column '{cid_column}' not found in obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    pairs = (
        adata.obs[[cluster_header, cid_column]]
        .astype(str)
        .drop_duplicates()
        .sort_values([cluster_header, cid_column])
        .reset_index(drop=True)
    )

    out = pd.DataFrame(
        {
            "cluster_name": pairs[cluster_header],
            "skos": default_skos,
            "manual_mapped_cid": "",
            "cell_ontology_id": pairs[cid_column],
        }
    )
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("h5ad", help="Path to input .h5ad file")
    p.add_argument("cluster_header", help="Name of cluster column in adata.obs")
    p.add_argument(
        "-o", "--output",
        default=None,
        help="Output CSV path (default: <h5ad-stem>_cluster_cid_mapping.csv)",
    )
    p.add_argument(
        "--cid-column",
        default="cell_type_ontology_term_id",
        help="obs column holding the current cell-ontology term ID "
             "(default: cell_type_ontology_term_id)",
    )
    p.add_argument(
        "--default-skos",
        default="",
        choices=["", "exact", "related"],
        help="Value to pre-populate in the skos column (default: blank)",
    )
    args = p.parse_args()

    out_path = (
        Path(args.output)
        if args.output
        else Path(args.h5ad).with_suffix("").with_name(
            Path(args.h5ad).stem + "_cluster_cid_mapping.csv"
        )
    )

    df = build_mapping(
        h5ad_path=args.h5ad,
        cluster_header=args.cluster_header,
        cid_column=args.cid_column,
        default_skos=args.default_skos,
    )
    df.to_csv(out_path, index=False)
    print(f"Wrote {len(df)} rows to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
