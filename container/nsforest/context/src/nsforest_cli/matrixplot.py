# container/nsforest/context/src/nsforest_cli/matrixplot.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import scanpy as sc


def _slice_leaves(
    leaves: list[str],
    *,
    leaf_range: Optional[str] = None,
    leaf_indices: Optional[Sequence[int]] = None,
) -> list[str]:
    if leaf_indices:
        n = len(leaves)
        idx = [i for i in leaf_indices if -n <= i < n]
        return [leaves[i] for i in idx]
    if leaf_range:
        try:
            start_s, end_s = leaf_range.split(":", 1)
            start = int(start_s) if start_s.strip() != "" else None
            end = int(end_s) if end_s.strip() != "" else None
            return leaves[slice(start, end)]
        except Exception:
            return leaves
    return leaves


def _resolve_symbol_map(adata, genes: List[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if hasattr(adata, "var") and "gene_symbols" in adata.var.columns:
        subset = [g for g in genes if g in adata.var.index]
        if subset:
            symbols = adata.var.loc[subset, "gene_symbols"].astype(str).to_dict()
            for g in subset:
                s = symbols.get(g)
                if s and s != "nan":
                    mapping[g] = s
    try:
        from nsforest_cli import ensembl_lookup as el
        for fn_name in ("lookup_symbols", "map_symbols", "map_ensembl_ids", "get_symbols"):
            if hasattr(el, fn_name):
                fn = getattr(el, fn_name)
                missing = [g for g in genes if g not in mapping]
                if missing:
                    try:
                        extra = fn(missing)
                        if isinstance(extra, dict):
                            mapping.update({k: v for k, v in extra.items() if v})
                    except Exception:
                        pass
                break
    except Exception:
        pass
    for g in genes:
        mapping.setdefault(g, g)
    return mapping


def matrixplot_run(
    *,
    h5ad_in: Path,
    results_csv: Path,
    label_key: str,
    png_out: Optional[Path] = None,
    svg_out: Optional[Path] = None,
    leaf_range: Optional[str] = None,
    leaf_indices: Optional[Sequence[int]] = None,
) -> None:
    """Matrix (dot) plot using nsforest.pl.matrixplot with log1p + symbols, options-only API.
    """
    import ast
    import numpy as np
    import pandas as pd
    from scipy import sparse as sp
    import nsforest as ns

    adata = sc.read_h5ad(str(h5ad_in))
    if label_key not in adata.obs:
        raise ValueError(f"obs['{label_key}'] not found in {h5ad_in}")
    adata.obs[label_key] = adata.obs[label_key].astype("category")

    dendro_key = f"dendrogram_{label_key}"
    if dendro_key not in adata.uns or "categories_ordered" not in adata.uns.get(dendro_key, {}):
        sc.tl.dendrogram(adata, groupby=label_key)
    leaves = list(adata.uns[dendro_key]["categories_ordered"])

    sliced = _slice_leaves(leaves, leaf_range=leaf_range, leaf_indices=leaf_indices)
    if clusters:
        want = set(clusters)
        dendrogram = [c for c in sliced if c in want]
    else:
        dendrogram = sliced
    if not dendrogram:
        raise ValueError("No clusters selected for plotting after subsetting.")

    adata = adata[adata.obs[label_key].isin(dendrogram)].copy()
    adata.obs[label_key] = adata.obs[label_key].cat.set_categories(dendrogram, ordered=True)

    df = pd.read_csv(results_csv)
    for col in ("cluster_header", "NSForest_markers"):
        if col not in df.columns:
            raise ValueError(f"Column '{col}' missing in {results_csv}")
    df = df[df["cluster_header"].isin(dendrogram)].copy()
    df[cluster_col] = df["cluster_header"].astype("category").cat.set_categories(dendrogram)
    df = df.sort_values("cluster_header")

    if top_n is not None and top_n > 0:
        def _trim(v: str) -> str:
            if not isinstance(v, str):
                return v
            s = v.strip()
            if s.startswith("[") and s.endswith("]"):
                try:
                    lst = [str(x).strip() for x in ast.literal_eval(s)]
                except Exception:
                    lst = [p.strip() for p in s.strip("[]").split(",")]
            else:
                lst = [p.strip() for p in s.split(";")]
            return ";".join(lst[:top_n])
        df["NSForest_markers"] = df["NSForest_markers"].map(_trim)

    def _split(val: str) -> list[str]:
        if not isinstance(val, str):
            return []
        s = val.strip()
        if s.startswith("[") and s.endswith("]"):
            try:
                lst = list(ast.literal_eval(s))
                return [str(x).strip() for x in lst if str(x).strip()]
            except Exception:
                pass
        return [g.strip() for g in s.split(";") if g.strip()]

    markers_dict: Dict[str, List[str]] = {
        row["cluster_header"]: _split(row["NSForest_markers"]) for _, row in df.iterrows() if row["cluster_header"] in dendrogram
    }

    union_order: List[str] = []
    seen = set()
    for cl in dendrogram:
        for g in markers_dict.get(cl, []):
            if g not in seen:
                seen.add(g)
                union_order.append(g)
    if not union_order:
        raise ValueError("No marker genes found to plot (empty union).")
    present = [g for g in union_order if g in adata.var_names]
    if not present:
        raise ValueError("None of the requested marker genes are present in adata.var_names.")

    sym_map = _resolve_symbol_map(adata, present)

    ad_for_plot = adata[:, present].copy()
    ad_for_plot.raw = None
    X = ad_for_plot.X
    if sp.issparse(X):
        ad_for_plot.X = sp.csr_matrix(np.log1p(X.toarray()))
    else:
        ad_for_plot.X = np.log1p(X)

    disp_names = [sym_map[g] for g in present]
    ad_for_plot.var_names = disp_names
    disp_set = set(disp_names)
    remapped_markers: Dict[str, List[str]] = {}
    for cl, genes in markers_dict.items():
        out: List[str] = []
        for g in genes:
            s = sym_map.get(g, g)
            if s in disp_set:
                out.append(s)
        remapped_markers[cl] = out

    import nsforest as ns
    ax = ns.pl.matrixplot(
        ad_for_plot,
        remapped_markers,
        label_key,
        dendrogram=dendrogram,
        save=False,
        output_folder=".",
        outputfilename_suffix=".",
    )
    fig = ax.get_figure() if hasattr(ax, "get_figure") else plt.gcf()

    if png_out:
        Path(png_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(png_out), bbox_inches="tight")
    if svg_out:
        Path(svg_out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(svg_out), bbox_inches="tight", format="svg")

    plt.close(fig)
    return None

