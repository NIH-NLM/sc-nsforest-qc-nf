# Why a consistent, containerized process — not "just NSForest"

## Thesis

The NSForest **algorithm** in this workflow is used **entirely unmodified**. That is a
deliberate design choice, not an oversight: the container `git clone`s the stock
https://github.com/NIH-NLM/NSForest and puts it on `PYTHONPATH` as-is
(`container/nsforest/Dockerfile`), and the code calls the library's own functions —
`ns.pp.prep_medians()`, `ns.pp.prep_binary_scores()`, `nsforesting.NSForest()`,
`ns.pp.dendrogram()`, `ns.pl.dotplot/stackedviolin/matrixplot/boxplot/scatter_w_clusterSize()`.

So the critique "this was all possible with NSForest without modification" is **true at the
algorithm level — and beside the point.** The value was never a change to the marker-selection
math. The value is everything required to run that algorithm **reproducibly, at scale, over
many heterogeneous real-world datasets, with QC-defined cohorts, comparable figures, and
silhouette-based cluster validation coupled in.** The stock `DEMO_NS-Forest_workflow` notebook
demonstrates that the algorithm *runs once, by hand, on one clean dataset*. It is not a process.

This document is the argument for everyone — production runs **and** individual experiments —
using the *same* containerized CLI, so that results are comparable, provenance-stamped, and
reproducible by construction.

---

## 1. What the NSForest package gives you vs. what this adds

| Concern | NSForest package alone (DEMO notebook) | sc-nsforest-qc-nf / `nsforest-cli` adds |
|---|---|---|
| Scope | One dataset, run by hand in a notebook | Many datasets from a `datasets_csv`, with a curation gate (`reference` = yes/no/unk vs exclude/merge/…) — `main.nf:57-68` |
| Execution | `NSForest()` called once, serially, over all clusters in one process | Scatter/gather: NSForest parallelized by cluster batch, partial CSVs merged — `main.nf:137-160`, `run_nsforest.py` |
| prep ↔ run | medians/binary computed in-memory, used immediately | Decoupled: prep computed once, serialized to CSV/pkl, re-attached to `varm` per batch — `prep.py`, `run_nsforest.py:47-89` |
| Cohort definition (QC) | none — feed it whatever cells you have | Ontology-ID filtering (UBERON tissue / MONDO+PATO disease / HsapDv age), per-row ID intersection, min-cluster-size, before/after dendrograms + stats — `filter_adata.py` |
| Reproducibility | depends on live env + network | Pinned container, gene map **baked at build** (no runtime network), provenance-stamped filenames `organ_author_journal_year_clusterheader_embedding_vid` — `common_utils.py:14-20` |
| Gene identifiers | uses whatever `var_names` are | Every artifact emitted in **both** ENSG and gene-symbol flavors from a controlled mapping — `merge_nsforest_results.py:33-56`, `gene_mapping_utils.py` |
| Cluster validation | none | **scsilhouette** scoring + QC viz coupled to the same run (see §3) |
| Ops | manual rerun | Nextflow `-resume`, per-process memory, trace/timeline/DAG, executor configs for Mac/Biowulf/AWS/CloudOS — `configs/` |

### Defensive fixes that make unattended runs possible

These only surface when NSForest is driven **unattended over messy CellxGene/Seurat data**. A
bare `NSForest()` call dies on each; the wrapper handles them so a 40-dataset run completes
without a human in the loop:

- **varm index doesn't persist across anndata versions** → DataFrames rebuilt with explicit ENSG
  index — `prep.py` (`_varm_to_df`).
- **varm alignment** → reindex to `var_names`, loud failure on real mismatch instead of silent
  NaN corruption — `run_nsforest.py:55-85`.
- **`ns.pp.dendrogram` hardcodes `use_rep="X_pca"`** → alias the dataset's real embedding into
  `X_pca`, three-level fallback for "negative distances" / low-dimensional embeddings —
  `filter_adata.py:74-145`.
- **scanpy 1.9.6 + matplotlib 3.8 `save="svg"` crash** → render via nsforest, save the figure
  ourselves; force `Agg` for headless containers — `filter_adata.py:15-18, 92-106`.
- **Dense `X` from large CXG datasets** → auto-convert to CSR (~10× memory) — `common_utils.py:28-43`.
- **Seurat-converted h5ads with positional `var.index`** → detect and promote the real gene-ID
  column — `filter_adata.py:387-409`.
- **anndata reserved `_index`** and **quotes/spaces/slashes in cluster names** (break filenames +
  CloudOS bash eval) → sanitized before write — `filter_adata.py:416-424, 542-553`.

---

## 2. The visualization layer — a consistent figure grammar

Plotting is not incidental here; it is a **standardized, comparable figure set** produced the same
way for every dataset and every person who uses the container. That comparability *is* the product.

Produced today by `nsforest-cli plots` (`plots.py`) and siblings, all gene-symbol-labeled and
dendrogram-ordered:

- **Boxplots** of `f_score`, `precision`, `recall`, `onTarget` — interactive **HTML** + publication
  **SVG** (`plots.py:58-64`).
- **Scatter vs. cluster size** for the same four metrics — HTML + SVG (`plots.py:66-72`).
- **Dotplot**, **stacked violin**, **matrix plot** of the NSForest-selected markers per cluster —
  each emitted in **raw** and **`standard_scale='var'` (scaled)** variants (`plots.py:83-105`).
- **Median / binary-score histograms** — `plot_histograms.py`.
- **Dendrogram** with auto-sized figure and the same negative-distance fallbacks — `dendrogram.py`.

Because every figure comes from one code path with provenance baked into the filename, a dotplot
from Anne's production run and a dotplot from a student's experiment are **directly comparable** —
same normalization, same gene-symbol mapping, same ordering, same file naming. That is exactly what
you cannot guarantee when everyone hand-rolls `sc.pl.dotplot` in their own notebook.

**This is the core argument for adoption:** the container is not "the production pipeline." It is the
*shared instrument*. Running an experiment through it costs nothing extra (every step is a
standalone `nsforest-cli` subcommand — see §4) and buys you results that line up with everyone
else's.

---

## 3. Why the scsilhouette + NSForest coupling matters

The workflow runs **scsilhouette** on the same filtered object that feeds NSForest and joins the two
in the same run — `compute_silhouette_process` plus `viz_summary`, `viz_distribution`,
`viz_2D_projection`, and `compute_summary_stats`, all joined against the merged NSForest results
(`main.nf:168-203`).

That coupling answers a question NSForest alone cannot: **are the clusters NSForest is finding
markers for actually well-separated?** A high `f_score` on a cluster with a poor silhouette is a
different claim than a high `f_score` on a cohesive cluster. Keeping the marker result and the
cluster-quality score in the same provenance-stamped run is what makes the marker set
*interpretable* rather than just *computed*. Experiments should carry this same coupling, not drop
it because they are "just exploring."

---

## 4. It already stands alone outside Nextflow

Nextflow is one caller of the container, not the thing itself. Every step is an independent
Typer subcommand (`main.py`):

```
nsforest-cli filter-adata …
nsforest-cli dendrogram …
nsforest-cli prep …               (medians + binary scores, one pass)
nsforest-cli run-nsforest …        nsforest-cli merge-nsforest-results …
nsforest-cli plots …               nsforest-cli plot-histograms …
nsforest-cli cluster-stats …       nsforest-cli cluster-cid-mapping …
nsforest-cli concat-h5ad …
```

So an experimenter who does **not** want the full production DAG can still `docker run` the same
image and call just the pieces they need — on their own h5ad, at the bench — and get outputs
that are byte-for-byte comparable in structure to production. Standing the CLI up as a
first-class, separately-documented tool (independent of this repo's `main.nf`) is a small
packaging step, not a rewrite.

---

## 5. Proposed next capability: focused, ad-hoc subtree visualization

**Motivation.** Today `plots` renders the full cluster set with NSForest-selected markers. For
exploration and for the ISMB-style story (endometriosis epithelial cells), we want to **zoom in**:
take one or two levels of the dendrogram, restrict to that subtree, and render dotplot / violin /
matrix on a *specific biomarker combination* — and to **compare two biomarker combinations** side
by side over the same subtree.

This is a new **input path** (user-supplied genes + a cluster subset), reusing the existing,
already-hardened plotting code. Proposed subcommands:

### 5a. `nsforest-cli dendrogram-subtree`
Cut the dendrogram at N levels and emit, per subtree: the member cluster list (CSV) and a focused
dendrogram SVG. Implementation reuses the linkage already computed in
`adata.uns['dendrogram_<header>']` (scipy `fcluster` on the stored linkage); no new algorithm.

```
nsforest-cli dendrogram-subtree \
  --h5ad adata_filtered_….h5ad --cluster-header <col> --levels 2 \
  --organ … --first-author … --journal … --year … --embedding … --dataset-version-id …
# → subtree_<k>_clusters_<prefix>.csv, dendrogram_subtree_<k>_<prefix>.svg
```

### 5b. `nsforest-cli focused-plots`
Render dotplot / stacked-violin / matrix restricted to a **cluster subset** (a subtree from 5a, or
an explicit `--clusters` list) and a **gene set** (an explicit biomarker combination, symbols or
ENSG). Reuses `ns.pl.dotplot / stackedviolin / matrixplot` exactly as `plots.py` does, so
normalization / scaling / gene-symbol labeling stay identical to production figures.

```
nsforest-cli focused-plots \
  --h5ad adata_filtered_….h5ad --cluster-header <col> \
  --clusters "epithelial_A,epithelial_B,epithelial_C"  \
  --genes "EPCAM,KRT8,KRT18,PGR,ESR1"                    \
  --label epithelial_combo1 …
# → dotplot_epithelial_combo1_<prefix>.svg (+ _scaled), stackedviolin_…, matrixplot_…
```

### 5c. Biomarker-combination comparison
Two gene combinations over the **same** subtree, emitted as a matched pair (and, optionally, a
combined figure) so combination A vs combination B is a like-for-like visual comparison — the
epithelial-cell comparison from the ISMB talk, generalized:

```
nsforest-cli focused-plots \
  --h5ad … --cluster-header <col> --clusters-from subtree_3_clusters_<prefix>.csv \
  --genes-a "EPCAM,KRT8,KRT18"  --label-a comboA \
  --genes-b "PGR,ESR1,WNT7A"    --label-b comboB
# → dotplot_comboA_<prefix>.svg, dotplot_comboB_<prefix>.svg, … (matched violin/matrix pairs)
```

**Why this belongs in the container and not a one-off notebook:** the comparison is only meaningful
if both combinations are rendered with identical normalization, ordering, and symbol mapping over
an identically-defined cell subset. That guarantee is exactly what the shared container provides and
what an ad-hoc script does not. Precompute these subtree/focused views ahead of time in production
so reviewers and experimenters browse a consistent gallery rather than regenerating figures by hand.

---

## 6. One-line rebuttal

> The NSForest algorithm is used unmodified — on purpose. What this adds is the *process* around it:
> ontology-driven QC to define cohorts, per-batch parallel execution, a pinned reproducible container
> with provenance-stamped dual gene-ID outputs, silhouette-based cluster validation coupled to every
> run, a standardized comparable figure set, and ~a dozen defensive fixes that let a bare `NSForest()`
> survive real CellxGene data. The DEMO notebook proves the algorithm runs once by hand; this proves
> it runs reproducibly, comparably, and unattended across many datasets — which is what makes everyone
> using the same instrument worthwhile, in experiments as much as in production.
