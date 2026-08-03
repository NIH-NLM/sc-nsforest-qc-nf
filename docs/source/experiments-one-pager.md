# Run your experiments on the shared NSForest-CLI — without touching production

**Why use it:** the container is a *shared instrument*, not just the production pipeline.
Running your experiment through it costs nothing extra and buys results that are directly
comparable to everyone else's — same environment, same gene-symbol mapping, same normalization,
same figure grammar, same provenance-stamped filenames. A dotplot you make lines up with a
dotplot from production because it came from the *same code path*. That comparability is the
whole point, and a hand-rolled notebook cannot guarantee it.

You also inherit, for free:
- **Reproducible environment** — one pinned container; no `conda`/`pip` drift to debug.
- **Baked-in gene mapping** — ENSG ↔ symbol works offline, identically for everyone.
- **QC coupled in** — silhouette scores + cluster validation travel with the marker results.
- **Provenance by construction** — every file is named `organ_author_journal_year_clusterheader_embedding_vid`.

---

## Quickstart — no Python setup, just Docker

Every pipeline step is a standalone subcommand. Mount your data, call the piece you need:

```bash
docker run --rm -v "$PWD":/data \
  ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest:latest \
  nsforest-cli --help
```

```bash
# Example: NSForest on your own filtered h5ad
docker run --rm -v "$PWD":/data \
  ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest:latest \
  nsforest-cli run-nsforest \
    --h5ad-path /data/adata_filtered_myexpt.h5ad \
    --medians-csv /data/medians_ensg_myexpt.csv \
    --binary-scores-csv /data/binary_scores_ensg_myexpt.csv \
    --cluster-header author_cell_type \
    --organ endometrium --first-author me --journal expt --year 2026
```

Used standalone like this, the CLI has **no publish path at all** — it cannot write to the
production graph. This is the safest way to experiment.

---

## The one rule: do not populate the production graph

The knowledge graph (`NIH-NLM/cell-kn`) is populated by exactly one process,
`publish_results_process`, and **it only runs when `--github_token` is set** (`main.nf:266`).
If you run the full Nextflow workflow, keep your experiment isolated with these knobs:

| Goal | What to do |
|---|---|
| **Never push to the graph** | **Omit `--github_token`.** Without it, the publish step is skipped and logs `WARNING: --github_token not set -- skipping publish step` (`main.nf:278`). Nothing reaches `NIH-NLM/cell-kn`. |
| **Keep local outputs separate** | Set `--outdir /your/experiment/dir` (default is `results/`). |
| **Publish, but to *your* space** | Override **both** `--publish_repo <you>/<your-fork>` and `--publish_dest_dir <your/path>`. Never leave them at the defaults (`NIH-NLM/cell-kn`, `data/prod/...`). |
| **Your own dataset list** | Point `--datasets_csv` at your CSV, not the production one. |
| **Never edit `main`** | **Fork the repo.** Make param/config/branch changes on your fork. Production runs from `main`; your fork can't affect it. |

### Safe experiment run (does not touch production)

```bash
nextflow run . \
  --datasets_csv   my_experiment.csv \
  --organ          endometrium \
  --uberon_json    data/uberon_endometrium.json \
  --disease_json   data/disease_endometriosis.json \
  --hsapdv_json    data/hsapdv_reproductive_age.json \
  --outdir         ./my_experiment_results
  # NOTE: no --github_token  → publish step is skipped, graph untouched
```

### If you fork and want your own publish target

```bash
nextflow run . \
  ... \
  --github_token     "$MY_TOKEN" \
  --publish_repo     me/my-cell-kn-fork \
  --publish_dest_dir data/experiments/me/endometrium/2026-06 \
  --outdir           ./my_experiment_results
```

---

## Checklist before you run

- [ ] I forked the repo (production runs from `main`; I work on my fork).
- [ ] I am **not** passing `--github_token`, **or** I redirected `--publish_repo` **and**
      `--publish_dest_dir` to my own space.
- [ ] `--outdir` points to my own directory.
- [ ] `--datasets_csv` is my list, not production's.
- [ ] I'm using the published image tag, so my environment matches everyone else's.

Follow this and your experiments are byte-for-byte comparable to production results while being
provably unable to alter the graph.
