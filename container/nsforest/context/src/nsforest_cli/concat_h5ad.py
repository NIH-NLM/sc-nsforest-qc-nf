"""
Concatenate multiple h5ad files into one, streaming on disk.

Apriori step before running the sc-nsforest-qc-nf workflow on a single
combined dataset. Uses anndata.experimental.concat_on_disk so RAM stays
flat regardless of input size.
"""

import time
from pathlib import Path

from anndata.experimental import concat_on_disk

from .common_utils import logger, log_section, setup_file_logging


def run_concat_h5ad(input_paths, output_path, join="outer", label="source_file"):
    """
    Concatenate h5ad files along the obs axis, streaming to disk.

    Args:
        input_paths: list of paths to input h5ad files
        output_path: path to write the concatenated h5ad
        join: 'outer' (union of vars, missing filled with zero) or 'inner'
        label: obs column name recording which input file each cell came from
    """
    setup_file_logging("concat_h5ad")
    log_section("NSForest: concat-h5ad (streaming)")

    if len(input_paths) < 2:
        raise ValueError(f"Need at least 2 input files; got {len(input_paths)}")

    in_paths = [Path(p) for p in input_paths]
    for p in in_paths:
        if not p.exists():
            raise FileNotFoundError(f"Input h5ad not found: {p}")
        if p.stat().st_size == 0:
            raise ValueError(f"Input h5ad is empty: {p}")

    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    keys = [p.stem for p in in_paths]
    if len(set(keys)) != len(keys):
        keys = [f"{i}_{p.stem}" for i, p in enumerate(in_paths)]

    logger.info(f"Inputs ({len(in_paths)}):")
    for k, p in zip(keys, in_paths):
        logger.info(f"  [{k}] {p}  ({p.stat().st_size / 1e9:.2f} GB)")
    logger.info(f"Output: {out_path}")
    logger.info(f"join={join}  label={label}")

    t0 = time.time()
    concat_on_disk(
        in_files=[str(p) for p in in_paths],
        out_file=str(out_path),
        axis=0,
        join=join,
        label=label,
        keys=keys,
        index_unique="-",
        merge="same",
        uns_merge="same",
    )
    elapsed = time.time() - t0
    logger.info(f"concat_on_disk complete in {elapsed:.1f}s")
    logger.info(f"Output size: {out_path.stat().st_size / 1e9:.2f} GB")
    
