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
    // val meta
    // files from all the files -- it all flat -- so I have to encode the context into the filename that is flattened
    tuple val(metas), file(all_files_ch)

    output:
    path "publish_report_${meta.organ}_${meta.first_author}_${meta.year}.txt", emit: report

    script:
    def today     = new java.text.SimpleDateFormat("yyyy-MMM-dd").format(new Date()).toLowerCase()
    def organSlug = metas.organ.replace('_', '-')
    def branch    = "${today}-${organSlug}-${metas.first_author}-${metas.year}-sc-nsforest-qc-nf"
    def label     = "outputs_${metas.organ}_${metas.first_author}_${metas.year}"
    def repo_url  = "https://\${GITHUB_TOKEN}@github.com/NIH-NLM/cell-kn.git"
    def report    = "publish_report_${metas.organ}_${metas.first_author}_${metas.year}.txt"
    """
    ls -lh 

    export GITHUB_TOKEN="${params.github_token}"

    echo "=========================================="
    echo " organ  : ${metas.organ}"
    echo " author : ${metas.first_author}"
    echo " year   : ${metas.year}"
    echo " branch : ${branch}"
    echo "=========================================="

#    git clone --depth 1 ${repo_url} cell-kn
#    cd cell-kn/prod/data/{metas.organ}/sc-nsforest-qc-nf
    
#    git config user.email "adeslatt@scitechcon.org"
#    git config user.name  "adeslatt"
#    git checkout -b ${branch}

#    git commit -m "workflow: publish ${metas.organ} ${metas.first_author} ${metas.year} results (${today})

#    git push origin ${branch}
    """
}
