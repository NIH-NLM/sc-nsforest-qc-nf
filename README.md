# sc-nsforest-qc-nf

[![Documentation Status](https://github.com/NIH-NLM/sc-nsforest-qc-nf/actions/workflows/docs.yml/badge.svg)](https://nih-nlm.github.io/sc-nsforest-qc-nf/)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Nextflow pipeline for NSForest marker gene discovery and silhouette score quality control of single-cell RNA-seq data.

`sc-nsforest-qc-nf` orchestrates parallel execution of [NSForest](https://github.com/JCVenterInstitute/NSForest) marker discovery and [scsilhouette](https://github.com/NIH-NLM/scsilhouette) clustering quality control across multiple datasets and organs, with ontology-based cell filtering driven by [cellxgene-harvester](https://github.com/NIH-NLM/cellxgene-harvester) outputs.

## Documentation

Full documentation is available at: **[https://nih-nlm.github.io/sc-nsforest-qc-nf/](https://nih-nlm.github.io/sc-nsforest-qc-nf/)**

## Features

- **Parallelized NSForest** marker gene discovery — scatter/gather by cluster for large-scale datasets
- **Silhouette QC** integrated with NSForest F-scores for comprehensive clustering assessment
- **Ontology-based filtering** — tissue (UBERON), disease (PATO:0000461), and age filters using the same logic as cellxgene-harvester
- **Modular architecture** — independent `modules/nsforest` and `modules/scsilhouette` process libraries
- **CloudOS-ready** — designed for scalable execution on cloud infrastructure
- **Standardized outputs** — consistent directory structure across all datasets and organs

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04.0
- [Docker](https://www.docker.com/) (containers pulled automatically)

## Containers

The pipeline uses two containers:

| Process group | Container |
|---|---|
| NSForest (all `modules/nsforest` processes) | `ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest:latest` |
| scsilhouette (all `modules/scsilhouette` processes) | `ghcr.io/nih-nlm/scsilhouette:1.0` |

Pull manually if needed:
```bash
docker pull ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest:latest
docker pull ghcr.io/nih-nlm/scsilhouette:1.0
```

## Installation

```bash
git clone https://github.com/NIH-NLM/sc-nsforest-qc-nf.git
cd sc-nsforest-qc-nf
```

No Python installation is required — all computation runs inside containers.

## Inputs

### Datasets CSV (`homo_sapiens_{organ}_harvester_final.csv`)

Produced by [cellxgene-harvester](https://github.com/NIH-NLM/cellxgene-harvester). One row per dataset. Required columns:

| Column | Description |
|---|---|
| `h5ad_file` | Path to input h5ad file (CellxGene format) |
| `first_author` | First author surname — used in output directory naming |
| `year` | Publication year — used in output directory naming |
| `author_cell_type` | `obs` column name for cluster labels |
| `embedding` | Embedding key (e.g. `X_umap`) |
| `disease` | Disease state label for annotation output |
| `filter_normal` | `True` / `False` — apply tissue + disease + age filtering |
| `reference` | Processing flag — see below |

**`reference` column values:**

| Value | Behaviour |
|---|---|
| `yes`, `no`, `unk` | Row is processed |
| `question`, `merge` | Row is **skipped** — logged at INFO level |

### UBERON JSON (`uberon_{organ}.json`)

Produced by `cellxgene-harvester resolve-uberon`. Encodes the full anatomical unit for the organ being processed — for example, `heart` includes heart, pericardium, and all UBERON descendant terms. The same file is used for every dataset in the run and drives:

1. **Tissue filter** — `tissue_ontology_term_id.isin(obo_ids)`
2. **Disease filter** — `disease_ontology_term_id == PATO:0000461`
3. **Age filter** — `development_stage` text parsed, cells with age ≥ `--min_age` retained

These filters are applied identically by both `filter-adata` (nsforest-cli) and `compute-silhouette` (scsilhouette).

### Input file location

Both input files are deposited manually after human review into:
```
cell-kn/data/prod/{organ}/cellxgene-harvester/
    homo_sapiens_{organ}_harvester_final.csv
    uberon_{organ}.json
```

## Quick Start

```bash
nextflow run main.nf \
    --datasets_csv homo_sapiens_kidney_harvester_final.csv \
    --organ kidney \
    --uberon_json uberon_kidney.json \
    --outdir results/kidney
```

### All parameters

| Parameter | Default | Description |
|---|---|---|
| `--datasets_csv` | required | Path to `homo_sapiens_{organ}_harvester_final.csv` |
| `--organ` | required | Organ label (e.g. `kidney`, `heart`) — used in output directory naming |
| `--uberon_json` | required | Path to `uberon_{organ}.json` |
| `--min_age` | `15` | Minimum donor age in years for adult cell filtering |
| `--min_cluster_size` | `5` | Minimum cells per cluster — smaller clusters dropped and logged |
| `--outdir` | `./results` | Output directory |
| `--publish_mode` | `copy` | Nextflow `publishDir` mode |

## Workflow Architecture

```
homo_sapiens_{organ}_harvester_final.csv
uberon_{organ}.json
        │
        ▼
[0] filter_adata         ← tissue (UBERON) + disease (PATO) + age + min cluster size
        │
        ├──────────────────────────────────────────────────┐
        ▼                                                  ▼
[1] dendrogram                                   [10] compute_silhouette
[2] cluster_stats                                [11] viz_summary
        │                                        [11] viz_distribution
        ▼ scatter by cluster                     [11] viz_dotplot
[4] prep_medians ×N                              [12] compute_summary_stats
        │
        ▼ gather
[5] merge_medians
        │
        ▼ scatter by cluster
[6] prep_binary_scores ×N
        │
        ▼ gather
[6] merge_binary_scores
        │
        ├── [7] plot_histograms
        │
        ▼ scatter by cluster
[8] run_nsforest ×N
        │
        ▼ gather
[8] merge_nsforest_results
        │
        ▼
[9] plots
```

The scatter/gather pattern (steps 4–8) parallelizes the computationally intensive median, binary score, and NSForest steps independently across every cluster in every dataset. A dataset with 50 clusters runs 50 parallel jobs at each scatter stage.

## Module Structure

```
modules/
├── nsforest/
│   ├── filter_adata.nf          ← tissue/disease/age/min-cluster filtering
│   ├── dendrogram.nf            ← cluster dendrogram + order
│   ├── cluster_stats.nf         ← cluster cell counts (drives scatter)
│   ├── prep_medians.nf          ← median expression per cluster (parallel)
│   ├── merge_medians.nf         ← gather partial medians
│   ├── prep_binary_scores.nf    ← binary scores per cluster (parallel)
│   ├── merge_binary_scores.nf   ← gather partial binary scores
│   ├── plot_histograms.nf       ← non-zero median/binary score histograms
│   ├── run_nsforest.nf          ← NSForest per cluster (parallel)
│   ├── merge_nsforest_results.nf← gather NSForest results
│   └── plots.nf                 ← boxplots, scatter, expression plots
└── scsilhouette/
    ├── compute_scsilhouette.nf  ← silhouette scores + cluster summary
    ├── viz_summary.nf           ← silhouette + F-score summary plot
    ├── viz_dotplot.nf           ← UMAP/embedding coloured by silhouette
    ├── viz_distribution.nf      ← cluster size vs silhouette distribution
    └── compute_summary_stats.nf ← dataset-level aggregate statistics
```

## NSForest CLI

The `modules/nsforest` processes call `nsforest-cli`, a workflow-internal command-line wrapper around the [NSForest](https://github.com/JCVenterInstitute/NSForest) library from the J. Craig Venter Institute. It is not published as a standalone package — it is bundled inside the `ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest` container and is specific to this workflow.

For the underlying NSForest algorithm, marker gene selection methodology, and citation information, refer to the **[NSForest repository](https://github.com/JCVenterInstitute/NSForest)**.

## Output Structure

All outputs follow the pattern: `results/{outputs_{organ}_{first_author}_{year}}/`

```
results/
└── outputs_kidney_Lake_2023/
    ├── adata_filtered.h5ad
    ├── subclass.full_cluster_sizes_before_filter.csv
    ├── subclass.full_cluster_sizes.csv
    ├── subclass.full_cluster_order.csv
    ├── subclass.full_summary_normal.csv
    ├── subclass.full_cluster_statistics.csv
    ├── subclass.full_medians.csv
    ├── subclass.full_binary_scores.csv
    ├── subclass.full_results.csv
    ├── subclass.full_silhouette_scores.csv
    ├── subclass.full_cluster_summary.csv
    ├── subclass.full_annotation.json
    ├── subclass.full_dataset_summary.csv
    ├── subclass.full_silhouette_fscore_summary.html
    ├── subclass.full_silhouette_fscore_summary.svg
    ├── subclass.full_dotplot_X_umap.html
    ├── subclass.full_dotplot_X_umap.svg
    ├── subclass.full_distribution_log10.html
    ├── subclass.full_distribution_log10.svg
    └── hist_nonzero_*.svg
```

## Using on Apple Silicon (M1/M2/M3)

Both containers are built for `linux/amd64`. On Apple Silicon Macs, add the `--platform` flag when pulling:
```bash
docker pull --platform linux/amd64 ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest:latest
docker pull --platform linux/amd64 ghcr.io/nih-nlm/scsilhouette:1.0
```

Nextflow will handle this automatically when running the pipeline with Docker on Apple Silicon.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License — see the LICENSE file for details.

## Acknowledgments

Developed by the National Library of Medicine (NLM), National Institutes of Health (NIH).

Part of the Cell Knowledge Network project.

## Contact

For questions or issues, please [open an issue](https://github.com/NIH-NLM/sc-nsforest-qc-nf/issues) on GitHub.

## Related Projects

- [NSForest](https://github.com/JCVenterInstitute/NSForest) — Marker gene discovery (J. Craig Venter Institute)
- [scsilhouette](https://github.com/NIH-NLM/scsilhouette) — Silhouette score QC package
- [cellxgene-harvester](https://github.com/NIH-NLM/cellxgene-harvester) — Single-cell data aggregation from CellxGene
- [cell-kn](https://github.com/NIH-NLM/cell-kn) — NIH NLM Cell Knowledge Network

## Citation

If you use sc-nsforest-qc-nf in your research, please cite:
```
[Citation information will be added upon publication]
```
