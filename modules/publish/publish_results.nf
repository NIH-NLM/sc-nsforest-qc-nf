/**
 * Publish Results to cell-kn GitHub Repository
 *
 * Fires once per dataset. Reads files from S3 results path derived
 * from workflow.workDir (works on CloudOS and other platforms).
 *
 * Branch naming:
 *   {YYYY}-{mon}-{DD}-{organ}-{first_author}-{year}-sc-nsforest-qc-nf
 *   e.g. 2026-mar-06-skin-of-body-Wiedemann-2023-sc-nsforest-qc-nf
 *
 * Platform notes:
 *   CloudOS: publishDir files land at s3://.../jobs/{id}/results/results/{label}/
 *            derived by replacing /work with /results/results in workflow.workDir
 */
process publish_results_process {
    tag "publish_results_${meta.organ}_${meta.first_author}_${meta.year}"
    label 'publish'

    input:
    tuple val(meta), file(all_files_ch)

    output:
    path "publish_report_${meta.organ}_${meta.first_author}_${meta.year}.txt", emit: report

    script:
    def today        = new java.text.SimpleDateFormat("yyyy-MMM-dd").format(new Date()).toLowerCase()
    def organ        = meta.organ
    def first_author = meta.first_author
    def year         = meta.year
    def organSlug    = organ.replace('_', '-')
    def branch       = "${today}-${organSlug}-${first_author}-${year}-sc-nsforest-qc-nf"
    def label        = "outputs_${organ}_${first_author}_${year}"
    def repo_url     = "https://\${GITHUB_TOKEN}@github.com/NIH-NLM/cell-kn.git"
    def report       = "publish_report_${organ}_${first_author}_${year}.txt"
    def run_id       = 123456789
    def dest_dir     = "prod/data/${organ}/sc-nsforest-qc-nf/${run_id}"
    """
    ls -lh

    export GITHUB_TOKEN="${params.github_token}"

    echo "=========================================="
    echo " organ  : ${organ}"
    echo " author : ${first_author}"
    echo " year   : ${year}"
    echo " branch : ${branch}"
    echo "=========================================="

    echo "publish complete: ${organ} ${first_author} ${year} branch: ${branch}" > ${report}
    git clone --depth 1 ${repo_url} cell-kn
    cd cell-kn

    git config user.email "adeslatt@scitechcon.org"
    git config user.name  "adeslatt"
    git checkout -b ${branch}

    mkdir -p ${dest_dir}
    mkdir -p results
    
    cp -L ../*.html ../*.log ../*.svg ../*.pkl ../*.json ../*.csv ${dest_dir}/ 2>/dev/null || true

    git add ${dest_dir}/
    git commit -m "workflow: publish ${organ} ${first_author} ${year} results (${today})"

    git push --force origin ${branch}
    """
}
