/**
 * Publish Results to cell-kn GitHub Repository
 *
 * Fires once per dataset. Reads files from S3 results path derived
 * from workflow.workDir (works on CloudOS and other platforms).
 *
 * Branch naming:
 *   {YYYY}-{mon}-{DD}-{organ}-{first_author}-{year}-sc_nsforest_qc_nf
 *   e.g. 2026-mar-06-skin-of-body-Wiedemann-2023-sc_nsforest_qc_nf
 *
 * Platform notes:
 *   CloudOS: publishDir files land at s3://.../jobs/{id}/results/results/{label}/
 *            derived by replacing /work with /results/results in workflow.workDir
 */
process publish_results_process {
    tag "publish_results_${meta.organ}_${meta.first_author}_${meta.year}"
    label 'publish'

    input:
    tuple val(meta), path('*')

    output:
    path "publish_report_${meta.organ}_${meta.first_author}_${meta.year}.txt", emit: report

    script:
    def today             = new java.text.SimpleDateFormat("yyyy-MMM-dd").format(new Date()).toLowerCase()
    def organ             = meta.organ
    def first_author      = meta.first_author
    def year              = meta.year
    def vid               = meta.dataset_version_id[-6..-1]
    def organSlug         = organ.replace('_', '-')
    def branch            = "${today}-${organSlug}-${first_author}-${year}-${vid}-sc_nsforest_qc_nf"
    def label             = "outputs_${organ}_${first_author}_${year}"
    def repo_url          = "https://\${GITHUB_TOKEN}@github.com/NIH-NLM/cell-kn.git"
    def report            = "publish_report_${organ}_${first_author}_${year}.txt"
    def run_id            = meta.session_id
    def sc_nsforest_qc_nf = "sc-nsforest-qc-nf"
    def dest_dir          = "data/prod/${sc_nsforest_qc_nf}/${organ}/${organ}-${first_author}-${year}-${vid}/${run_id}/results"
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

    cp -L ../*.html ../*.log ../*.svg ../*.pkl ../*.json ../*.csv ${dest_dir}/ 2>/dev/null || true

    # Compress files larger than 50MB
    find ${dest_dir} -type f -size +50M | while read f; do
        tar czf "\${f}.tar.gz" -C "\$(dirname "\$f")" "\$(basename "\$f")"
        rm "\$f"
    done

    git add ${dest_dir}/
    git commit -m "workflow: publish ${organ} ${first_author} ${year} ${vid} results (${today})"

    git push --force origin ${branch}
    """
}
