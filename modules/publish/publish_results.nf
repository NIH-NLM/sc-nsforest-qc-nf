/**
 * Publish Results to cell-kn GitHub Repository
 *
 * Fires once per dataset. Reads files from S3 results path derived
 * from workflow.workDir (works on CloudOS and other platforms).
 *
 * Branch naming:
 *   {YYYY}-{mon}-{DD}}-{session id [-6..-1]}-${organSlug}-sc_nsforest_qc_nf
 *   one branch per organ 
 *   e.g. 2026-mar-06-abcde-kidney-sc_nsforest_qc_nf
 *
 * Platform notes:
 *   CloudOS: publishDir files land at s3://.../jobs/{id}/results/{label}/
 *            derived by replacing /work with /results in workflow.workDir
 */
process publish_results_process {
    tag "publish_results"
    label 'publish'

    input:
    tuple val(meta), path('*')

    output:
    path "publish_report.txt", emit: report

    script:
    def today     = new java.text.SimpleDateFormat("yyyy-MMM-dd").format(new Date()).toLowerCase()
    def organ     = params.organ
    def organSlug = organ.replace('_', '-')
    def firstAuth = meta.first_author
    def authSlug  = firstAuth.replace(' ', '-').replace(',', '')
    def year      = meta.year
    def journal   = meta.journal
    def vid       = meta.dataset_version_id.toString()[-6..-1]
    def vidFull   = meta.dataset_version_id.toString()   // the whole sha256 of the unfiltered h5ad
    def sid       = workflow.sessionId.toString()[-6..-1]
    def branch    = "${today}-${sid}-${organSlug}-sc_nsforest_qc_nf"
    def dest_dir  = params.publish_dest_dir ?: "data/prod/${organ}/sc-nsforest-qc-nf/results/${organ}-${firstAuth}-${journal}-${year}-${vidFull}"
    def repo      = params.publish_repo ?: 'NIH-NLM/cell-kn'
    def repo_url  = "https://\${GITHUB_TOKEN}@github.com/${repo}.git"
    """
    ls -lh

    export GITHUB_TOKEN="${params.github_token}"

    echo "=========================================="
    echo " organ  : ${organ}"
    echo " author : ${firstAuth}"
    echo " year   : ${year}"
    echo " vid    : ${vid}"
    echo " sid    : ${sid}"
    echo " branch : ${branch}"
    echo " dest   : ${dest_dir}"
    echo "=========================================="

    echo "publish complete: ${organ} ${firstAuth} ${year} ${vid} branch: ${branch}" > publish_report.txt

    git clone --depth 1 ${repo_url} publish-repo
    cd publish-repo

    git config user.email "sc-nsforest-qc-nf@noreply.github.com"
    git config user.name  "sc-nsforest-qc-nf"
    git fetch --depth=1 origin ${branch}:${branch} 2>/dev/null || true
    git checkout ${branch} 2>/dev/null || git checkout -b ${branch}

    mkdir -p "${dest_dir}"
    # the pinned tools end their file names with the last 6 characters of the id: the copies in GitHub carry all of it
    for f in ../*.html ../*.log ../*.svg ../*.png ../*.pkl ../*.json ../*.csv; do
        [ -e "\$f" ] || continue
        name=\$(basename "\$f" | sed -e "s/_${vid}\\./_${vidFull}./" -e "s/_${vid}_/_${vidFull}_/")
        cp -L "\$f" "${dest_dir}/\$name"
    done

    find "${dest_dir}" -type f -size +50M | while read f; do
        tar czf "\${f}.tar.gz" -C "\$(dirname "\$f")" "\$(basename "\$f")"
        rm "\$f"
    done

    # GitHub refuses a file over 100 MB: a file still above 99 MB is not pushed (it stays in the S3 results folder and in the manifest)
    find "${dest_dir}" -type f -size +99M | while read f; do
        echo "not pushed to GitHub, over 99 MB: \$(basename "\$f")" | tee -a ../publish_report.txt
        rm "\$f"
    done

    git add "${dest_dir}/"
    git commit -m "workflow: publish ${organ} ${firstAuth} ${year} ${vid} results (${today})"
    git push origin ${branch}
    """
}
