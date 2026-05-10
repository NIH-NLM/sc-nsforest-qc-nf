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
    --disease disease_normal.json \
    --hsapdv hsapdv_15.json
    --github_token {{$secret_github_token}}
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
docker pull --platform linux/amd64 ghcr.io/nih-nlm/scsilhouette:latest
```

Nextflow will handle this automatically when running the pipeline with Docker on Apple Silicon.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

The float32 calculation

3.2M cells × 20K genes = 3,200,000 × 20,000 = 6.4 × 10¹⁰ entries
float32 = 4 bytes per entry
6.4 × 10¹⁰ × 4 = 2.56 × 10¹¹ bytes = 256 GB
That's the dense matrix size — every cell × gene cell stored. The actual h5ad on disk is much smaller because scRNA matrices are >90% zeros and stored sparse (CSR/CSC: just the nonzeros + indices). Sparse in-memory is typically 10-25 GB for a 3M-cell dataset. The OOM happens when something quietly densifies a slice (e.g. .toarray(), certain scanpy functions, casting through pandas).

Quick reckoner: dense bytes = cells × genes × dtype_bytes. float64 doubles it; int32 is the same; bool is 1 byte.

2. Significant digits for F-β in NSForest

There are two different things people mean by "significant digits" here:

A. Numerical precision — basically irrelevant. F-β is computed from integer TP/FP/FN counts:

F_β = (1 + β²) · TP / ((1 + β²) · TP + β² · FN + FP)
NSForest uses β = 0.5 (precision-weighted). The only floating-point step is the final division, so float32 gives ~7 decimal digits of numerical accuracy. Not your concern.

B. Statistical precision — this is what actually matters and depends on cluster size N. Treating precision P and recall R as binomial proportions, the standard error of F-β by error propagation is:

| Symbol | Name | Meaning |
| ------ | ---- | ------- |
| σ | sigma (lowercase) | standard deviation of one variable |
| σ² | sigma squared | variance |
| Σ | sigma (uppercase) | summation operator |
| μ | mu | mean or expected value |
| β | beta | the parameter in F-β (user sets this; NS-Forest uses 0.5 |
| ∂ | partial | partial derivative |


σ²(F_β) ≈ (∂F/∂P)² · σ²(P) + (∂F/∂R)² · σ²(R)
where:
  σ²(P) = P(1−P) / N_pred_pos
  σ²(R) = R(1−R) / N_actual_pos
  ∂F/∂P = (1+β²) · R² / (β²·P + R)²
  ∂F/∂R = (1+β²) · β² · P² / (β²·P + R)²

Rule of thumb: meaningful digits ≈ ½ · log₁₀(N_smallest_cluster).

100 cells → ~1 sig digit (e.g., 0.7)
10,000 cells → ~2 sig digits (e.g., 0.74)
1M cells → ~3 sig digits (0.743)
So if your smallest cluster has 50-100 cells, reporting F-β to 4 decimals is theater.

The practical answer: don't derive it analytically. Bootstrap. Run NSForest's evaluation step on B=100 resamples (with replacement) of the test cells per cluster and report the 2.5/97.5 percentiles of F-β. That captures the actual uncertainty including any non-binomial structure (correlated genes, etc.) that the closed-form ignores. For NSForest specifically, the bootstrap is one wrapper around their existing per-cluster evaluation — much less work than stepping through the variance algebra and gives you proper CIs.

## AWS Instance-type suffix details

| Suffix | What changes | When to pick it |
| ————-- | ——————————–- | ——————————————- |
| r5 | Intel Xeon Skylake/cascade Lake, EBS-only, standard networking. Baseline | Default; what most docs assume |
| r5a | AMD EPYC instead of Intel.  ~10% cheaper.  Slightly lower single-core performance on some workloads | Cost savings when you don't depend on Intel specific instructions |
| r5ad | AMD + NVMe instance-store SSD attached locally. | AMD pricing + you want fast ephermeral scratch (data lost on stop) |
| r5b | Intel + EBS-optimized with much higher EBS bandwidth (up to 60 Gbps vs ~14 Gpbs on r5 | EBS throughput is your bottleneck == e.g., reading huge h5ad files off gp3 fast. |
| r5d | Intel + NVMe instance-store SSD | Need fast local scratch, prefer Intel |
| r5n | Intel + enhanced networking (up to 100 Gbps). | Lots of S3 / cross AZ traffic |
| r5dn | Intel + NVMe + enhanced networking | Both: heavy network and heavy local-scratch I/O |

## TMI and TLDR

Cross-AZ traffic — "AZ" stands for Availability Zone. An AWS region (e.g., us-east-1, "Northern Virginia") is split into several physically separate datacenters called Availability Zones (us-east-1a, us-east-1b, us-east-1c, ...). Each zone has its own power, cooling, and networking. Cross-AZ traffic means data moving between two AZs inside the same region — e.g., from your EC2 instance in us-east-1a to a database in us-east-1b. AWS charges $0.01 per GB in each direction for this, and latency is ~1-2 ms vs ~0.1 ms within a single AZ. By contrast, EC2 → S3 inside the same region is free and stays inside AWS's backbone. Practical implication: keep your interactive instance, your S3 bucket, and any other working storage in the same region, and ideally launch the instance in the AZ closest to where your bucket replicates.

AVX-512 — Advanced Vector eXtensions, 512-bit edition. A set of CPU instructions that Intel added to its Xeon server chips starting around 2017. The idea: instead of doing one floating-point multiply per cycle, the CPU does 16 multiplies at once on a 512-bit-wide register (16 × 32-bit floats = 512 bits). This is a form of SIMD (Single Instruction, Multiple Data). Numerical libraries like numpy, MKL, OpenBLAS can use AVX-512 to accelerate matrix multiplication, FFTs, and similar work — typically 1.5-2× faster than the previous AVX2 (256-bit). Why it matters for instance choice: AMD EPYC processors (the "a" suffix on AWS instances) lacked AVX-512 until 2022. So picking r5a over r5 saved ~10% on cost but gave up AVX-512 speedups for numerically intensive code. As of the latest AMD generation, this is no longer a meaningful gap — but legacy code paths and older AMD instance types still don't have it.

AMD EPYC — AMD's server CPU brand, the competitor to Intel's Xeon. "EPYC" is just the marketing name (read "epic"). Pronounced like the English word. EPYC chips are typically denser (more CPU cores per socket) and cheaper per core than equivalent Xeons, with strong memory bandwidth. The trade-off historically was per-core performance (Intel was a bit faster on single-threaded code) and instruction-set features (no AVX-512 until recently). On AWS, the a in r5a / m5a / c5a means "AMD EPYC inside instead of Intel Xeon" — same vCPU count, same RAM, ~10% cheaper. For most bioinformatics work (which is heavily parallel and bottlenecked by I/O or memory bandwidth, not single-thread speed), EPYC is fine.

NVMe — Non-Volatile Memory Express. A protocol for talking to SSDs directly over the PCIe bus (the same bus that connects GPUs and high-speed network cards), bypassing the older SATA disk interface. The result: NVMe SSDs deliver millions of operations per second and 3-7 GB/s of sequential throughput, vs SATA SSDs at ~500 MB/s. On AWS, "instance store" SSDs on d-suffixed instances (e.g., r5d, r5ad, m5d) are NVMe drives physically attached to the host machine, not over the network. Pros: very fast (great for scratch / temp / shuffle data, HDF5 random reads). Cons: ephemeral — when you stop or terminate the instance, the NVMe disk is wiped. So they're perfect for compute-and-discard workloads (like concat-h5ad on disposable infrastructure) and wrong for anything you need to persist. EBS, by contrast, is network-attached storage — slower per IOP but durable across instance lifecycle.

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
