"""
Concatenate multiple h5ad files into one, streaming on disk.

Apriori step before running the sc-nsforest-qc-nf workflow on a single
combined dataset. Uses anndata.experimental.concat_on_disk so RAM stays
flat regardless of input size. Inputs are pre-sanitized to move 1D entries
mistakenly stored in varm/obsm into var/obs, since concat_on_disk rejects
non-2D values in those mappings.
"""

import tempfile
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from anndata.experimental import concat_on_disk

from .common_utils import logger, log_section, setup_file_logging

def _is_1d(value) -> bool:
    """True if value is a 1D Series or 1D ndarray-like."""
    if isinstance(value, pd.Series):
        return True
    shape = getattr(value, "shape", None)
    return shape is not None and len(shape) == 1


def _coerce_to_series(value, target_index, name):
    """Return a Series aligned to target_index, or None if unfixable."""
    n = len(target_index)
    if isinstance(value, pd.Series):
        if len(value) != n:
            return None
        if value.index.equals(target_index):
            return value
        if set(value.index).issubset(set(target_index)) or set(target_index).issubset(set(value.index)):
            return value.reindex(target_index)
        return pd.Series(value.values, index=target_index, name=name)
    arr = np.asarray(value)
    if arr.ndim == 1 and len(arr) == n:
        return pd.Series(arr, index=target_index, name=name)
    return None


def _repair_axis_mapping(adata, mapping_attr: str, frame_attr: str, file_label: str) -> int:
    """
    Move 1D entries out of varm/obsm into var/obs. Returns count of repairs.

    Rules:
      - 1D entry, fits target axis -> move into frame
      - name already in frame with equal data -> drop misplaced copy
      - name already in frame with different data -> save as f"{name}_from_{mapping_attr}"
      - 1D entry, length mismatch / uncoercible -> drop with WARN
    """
    mapping = getattr(adata, mapping_attr)
    frame = getattr(adata, frame_attr)
    target_index = frame.index
    bad_keys = [k for k, v in mapping.items() if _is_1d(v)]
    n_changes = 0

    for key in bad_keys:
        value = mapping[key]
        try:
            series = _coerce_to_series(value, target_index, key)
        except Exception as e:
            logger.warning(
                f"[sanitize:{file_label}] {mapping_attr}['{key}']: coercion failed ({e!r}); dropping."
            )
            del mapping[key]
            n_changes += 1
            continue

        if series is None:
            logger.warning(
                f"[sanitize:{file_label}] {mapping_attr}['{key}']: length mismatch with "
                f"{frame_attr} (len={len(target_index)}); dropping."
            )
            del mapping[key]
            n_changes += 1
            continue

        if key in frame.columns:
            try:
                same = frame[key].equals(series.astype(frame[key].dtype, errors="ignore"))
            except Exception:
                same = False
            if same:
                logger.info(
                    f"[sanitize:{file_label}] {mapping_attr}['{key}']: already in "
                    f"{frame_attr}; dropping copy."
                )
                del mapping[key]
            else:
                new_key = f"{key}_from_{mapping_attr}"
                logger.warning(
                    f"[sanitize:{file_label}] {mapping_attr}['{key}']: conflicts with "
                    f"{frame_attr}['{key}']; moved to {frame_attr}['{new_key}']."
                )
                frame[new_key] = series
                del mapping[key]
        else:
            logger.info(
                f"[sanitize:{file_label}] {mapping_attr}['{key}']: 1D — moved into "
                f"{frame_attr}['{key}']."
            )
            frame[key] = series
            del mapping[key]

        n_changes += 1

    return n_changes


def _sanitize_h5ad(in_path: Path, out_path: Path) -> bool:
    """
    Load in_path, repair 1D entries in varm/obsm, write to out_path if changed.

    Returns True if any change was made (out_path was written).
    Returns False if no repairs needed (caller should use in_path directly).
    """
    file_label = in_path.name
    adata = ad.read_h5ad(in_path)
    n = 0
    n += _repair_axis_mapping(adata, "varm", "var", file_label)
    n += _repair_axis_mapping(adata, "obsm", "obs", file_label)
    if n == 0:
        logger.info(f"[sanitize:{file_label}] clean — no repairs.")
        return False
    logger.info(f"[sanitize:{file_label}] {n} repair(s); writing sanitized copy.")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)
    return True

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
    with tempfile.TemporaryDirectory(prefix="nsforest_concat_sanitize_") as td:
        td_path = Path(td)
        effective_paths: list[Path] = []
        for i, p in enumerate(in_paths):
            sanitized = td_path / f"{i:03d}_{p.stem}.h5ad"
            try:
                changed = _sanitize_h5ad(p, sanitized)
                effective_paths.append(sanitized if changed else p)
            except Exception as e:
                logger.error(
                    f"[sanitize:{p.name}] failed ({e!r}); falling back to original. "
                    f"concat_on_disk may surface the underlying error."
                )
                effective_paths.append(p)

        n_repaired = sum(1 for ep, op in zip(effective_paths, in_paths) if ep != op)
        logger.info(f"Sanitization summary: {n_repaired}/{len(in_paths)} files required repairs.")

        concat_on_disk(
            in_files=[str(p) for p in effective_paths],
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
    
