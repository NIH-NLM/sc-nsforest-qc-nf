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
    def today     = new java.text.SimpleDateFormat("yyyy-MMM-dd").format(new Date()).toLowerCase()
    def organSlug = meta.organ.replace('_', '-')
    def branch    = "${today}-${organSlug}-${meta.first_author}-${meta.year}-sc-nsforest-qc-nf"
    def label     = "outputs_${meta.organ}_${meta.first_author}_${meta.year}"
    def repo_url  = "https://\${GITHUB_TOKEN}@github.com/NIH-NLM/cell-kn.git"
    def report    = "publish_report_${meta.organ}_${meta.first_author}_${meta.year}.txt"
    """
    ls -lh 

    export GITHUB_TOKEN="${params.github_token}"

    echo "=========================================="
    echo " organ  : ${meta.organ}"
    echo " author : ${meta.first_author}"
    echo " year   : ${meta.year}"
    echo " branch : ${branch}"
    echo "=========================================="

    echo "publish complete: ${meta.organ} ${meta.first_author} ${meta.year} branch: ${branch}" > ${report}
    git clone --depth 1 ${repo_url} cell-kn
    cd cell-kn/prod/data/${meta.organ}/sc-nsforest-qc-nf
   
    git config user.email "adeslatt@scitechcon.org"
    git config user.name  "adeslatt"
    git checkout -b ${branch}

    git commit -m "workflow: publish ${meta.organ} ${meta.first_author} ${meta.year} results (${today})"

    git push origin ${branch}
    """
}
