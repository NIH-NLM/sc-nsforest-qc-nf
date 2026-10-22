#!/usr/bin/env python3
"""
Make a small h5ad with the real structure of a large one, for rehearsing the whole workflow locally.

Takes up to --per-cluster cells from every cluster, so every cluster (and so every subclass) is present, and keeps
all the genes, the obsm, the raw counts and the way the strings are stored. The file is read in backed mode, so the
large expression matrix is not loaded into memory.

  python scripts/make_rehearsal_subsample.py adata_sampled_AIBS_human_neocortex_plainstrings.h5ad \
      rehearsal_subsample.h5ad --cluster-header CrossArea_cluster --per-cluster 20

Written to be read with the same anndata as the input (and the image's anndata 0.9.2 reads it if the input was
saved with plain strings).
"""
import argparse

import anndata as ad
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("h5ad_in")
    parser.add_argument("h5ad_out")
    parser.add_argument("--cluster-header", default="CrossArea_cluster")
    parser.add_argument("--per-cluster", type=int, default=20, help="most cells taken from each cluster")
    parser.add_argument("--seed", type=int, default=10)
    args = parser.parse_args()

    backed = ad.read_h5ad(args.h5ad_in, backed="r")
    labels = backed.obs[args.cluster_header].astype(str).to_numpy()
    rng = np.random.default_rng(args.seed)

    keep = []
    for cluster in np.unique(labels):
        idx = np.flatnonzero(labels == cluster)
        keep.append(idx if len(idx) <= args.per_cluster else rng.choice(idx, args.per_cluster, replace=False))
    keep = np.sort(np.concatenate(keep))

    small = backed[keep].to_memory()
    backed.file.close()
    small.write_h5ad(args.h5ad_out)

    print(f"{small.n_obs:,} cells x {small.n_vars:,} genes, {len(np.unique(labels))} clusters -> {args.h5ad_out}")
    for header in ("CrossArea_subclass", "CrossArea_cluster"):
        if header in small.obs.columns:
            print(f"  {header}: {small.obs[header].nunique()} groups, smallest has {small.obs[header].value_counts().min()} cells")


if __name__ == "__main__":
    main()
