/**
 * Publish Results to cell-kn GitHub Repository
 *
 * Fires ONCE after ALL datasets in the CSV have been processed (or failed).
 * Uses .collect() sentinel signals from the terminal processes of both the
 * NSForest branch (plots_process) and the scsilhouette branch
 * (compute_summary_stats_process) to guarantee nothing runs until every
 * parallel job has finished or errored out.
 *
 * The process works directly from params.outdir — no re-staging of files.
 * publishDir has already written everything there.
 *
 * Published directory structure in cell-kn:
 * ------------------------------------------
 *   data/prod/{organ}/nsforest/{organ}_{first_author}_{year}/
 *     All NSForest outputs: filter stats, dendrogram, medians, binary scores,
 *     histograms, marker results, plots.
 *
 *   data/prod/{organ}/scsilhouette/{organ}_{first_author}_{year}/
 *     All scsilhouette outputs: scores CSV, cluster summary, annotation JSON,
 *     viz plots, dataset summary.
 *
 * Datasets are identified by signature files inside each outputs_* directory:
 *   NSForest     → presence of *_results.csv   (merge_nsforest_results output)
 *   scsilhouette → presence of *_silhouette_scores.csv (compute_silhouette output)
 * A single outputs_* directory will contain BOTH if both branches completed.
 *
 * Branch naming:
 *   workflow/{organ}_{YYYY-MM-DD}_{workflow.runName}
 *   e.g. workflow/kidney_2025-01-15_hungry_lovelace
 *
 * Required params:
 *   params.github_token   GitHub PAT with repo write access.
 *                         Pass via --github_token or a nextflow.config params block.
 *                         NEVER hardcode.
 *   params.outdir         Local results directory written by all publishDir steps.
 *
 * Input:
 *   organ              val — organ string (e.g. kidney)
 *   nsforest_signals   val — collected list of dataset labels from plots_process.
 *                      Used only as a completion gate — values are not read.
 *   silhouette_signals val — collected list of dataset labels from
 *                      compute_summary_stats_process.
 *                      Used only as a completion gate — values are not read.
 *
 * Output:
 *   publish_report.txt  Written to params.outdir. Records branch, PR URL,
 *                       and counts of published datasets.
 */
process publish_results_process {
    tag "publish_${organ}"
    label 'publish'

    publishDir "${params.outdir}",
        mode: params.publish_mode,
        pattern: "publish_report.txt"

    input:
    val  organ
    val  nsforest_signals
    val  silhouette_signals

    output:
    path "publish_report.txt", emit: report

    script:
    def today    = new java.text.SimpleDateFormat("yyyy-MM-dd").format(new Date())
    def branch   = "workflow/${organ}_${today}_${workflow.runName}"
    def outdir   = params.outdir
    def repo_url = "https://\${GITHUB_TOKEN}@github.com/NIH-NLM/cell-kn.git"
    """
    set -euo pipefail

    export GITHUB_TOKEN="${params.github_token}"
    export GH_TOKEN="${params.github_token}"

    echo "=========================================="
    echo " organ   : ${organ}"
    echo " branch  : ${branch}"
    echo " outdir  : ${outdir}"
    echo "=========================================="

    # ------------------------------------------------------------------
    # Clone cell-kn (shallow — only HEAD of main needed)
    # ------------------------------------------------------------------
    git clone --depth 1 ${repo_url} cell-kn
    cd cell-kn

    git config user.email "nextflow-publish@nih.gov"
    git config user.name  "Nextflow Publish Bot"
    git checkout -b ${branch}

    # ------------------------------------------------------------------
    # NSForest outputs
    # Identified by *_results.csv produced by merge_nsforest_results
    # Destination: data/prod/{organ}/nsforest/{organ}_{author}_{year}/
    # ------------------------------------------------------------------
    echo ""
    echo "--- NSForest ---"
    n_nsforest=0
    for src_dir in "${outdir}"/outputs_${organ}_*/; do
        [ -d "\$src_dir" ] || continue
        ls "\$src_dir"/*_results.csv 2>/dev/null || continue
        dataset_label=\$(basename "\$src_dir" | sed 's/^outputs_//')
        dest="data/prod/${organ}/nsforest/\${dataset_label}"
        mkdir -p "\$dest"
        cp -rp "\$src_dir"/* "\$dest/"
        echo "  + \$dataset_label"
        n_nsforest=\$(( n_nsforest + 1 ))
    done
    echo "  total: \$n_nsforest"

    # ------------------------------------------------------------------
    # scsilhouette outputs
    # Identified by *_silhouette_scores.csv from compute_silhouette
    # Destination: data/prod/{organ}/scsilhouette/{organ}_{author}_{year}/
    # ------------------------------------------------------------------
    echo ""
    echo "--- scsilhouette ---"
    n_silhouette=0
    for src_dir in "${outdir}"/outputs_${organ}_*/; do
        [ -d "\$src_dir" ] || continue
        ls "\$src_dir"/*_silhouette_scores.csv 2>/dev/null || continue
        dataset_label=\$(basename "\$src_dir" | sed 's/^outputs_//')
        dest="data/prod/${organ}/scsilhouette/\${dataset_label}"
        mkdir -p "\$dest"
        cp -rp "\$src_dir"/* "\$dest/"
        echo "  + \$dataset_label"
        n_silhouette=\$(( n_silhouette + 1 ))
    done
    echo "  total: \$n_silhouette"

    # ------------------------------------------------------------------
    # Stage and commit — idempotent if nothing changed
    # ------------------------------------------------------------------
    git add data/prod/${organ}/

    if git diff --cached --quiet; then
        echo ""
        echo "WARNING: nothing new to commit — outputs already up to date"
        cat > publish_report.txt << REPORT
organ:         ${organ}
branch:        ${branch}
run_name:      ${workflow.runName}
date:          ${today}
nsforest:      \${n_nsforest} datasets
scsilhouette:  \${n_silhouette} datasets
status:        nothing to commit — already up to date
pr_url:        (none opened)
REPORT
        exit 0
    fi

    git commit -m "workflow: publish ${organ} results (${today})

Organ         : ${organ}
Run           : ${workflow.runName}
Date          : ${today}
NSForest      : \${n_nsforest} datasets
scsilhouette  : \${n_silhouette} datasets

Auto-generated by sc-nsforest-qc-nf."

    # ------------------------------------------------------------------
    # Push and open PR
    # ------------------------------------------------------------------
    git push origin ${branch}

    gh pr create \\
        --repo NIH-NLM/cell-kn \\
        --base main \\
        --head "${branch}" \\
        --title "workflow: publish ${organ} results (${today})" \\
        --body "## Automated results publish

| Field | Value |
|---|---|
| **Organ** | ${organ} |
| **Run name** | ${workflow.runName} |
| **Date** | ${today} |
| **NSForest datasets** | \${n_nsforest} |
| **scsilhouette datasets** | \${n_silhouette} |

### Directory structure
\`\`\`
data/prod/${organ}/
├── nsforest/
│   └── {organ}_{first_author}_{year}/
└── scsilhouette/
    └── {organ}_{first_author}_{year}/
\`\`\`

### Review checklist
- [ ] Dataset count matches CSV input (minus skipped / failed rows)
- [ ] Spot-check silhouette scores and NSForest marker genes
- [ ] Confirm failed datasets (if any) are acceptable before merge
- [ ] cellxgene-harvester outputs and human annotation deposited separately

_Auto-opened by sc-nsforest-qc-nf._"

    PR_URL=\$(gh pr view --repo NIH-NLM/cell-kn "${branch}" --json url -q '.url' 2>/dev/null || echo "unknown")

    echo ""
    echo "PR opened: \$PR_URL"

    cat > publish_report.txt << REPORT
organ:         ${organ}
branch:        ${branch}
run_name:      ${workflow.runName}
date:          ${today}
nsforest:      \${n_nsforest} datasets
scsilhouette:  \${n_silhouette} datasets
status:        published
pr_url:        \$PR_URL
REPORT
    cat publish_report.txt
    """
}
