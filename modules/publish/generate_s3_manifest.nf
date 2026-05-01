/**
 * Generate S3 Manifest
 *
 * Runs once across all datasets. Collects all published output files
 * and writes master_s3_manifest.csv listing each filename and its
 * permanent S3 path in the CloudOS results bucket.
 *
 * Input:
 * ------
 * @param files:    All output files collected across all datasets
 * @param s3_base:  S3 results base path (workflow.workDir.parent/results)
 *
 * Output:
 * -------
 * @emit manifest: master_s3_manifest.csv
 */
process generate_s3_manifest_process {
    tag "generate_s3_manifest"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    path('*')
    val s3_base

    output:
    path "master_s3_manifest.csv", emit: manifest

    script:
    """
    nsforest-cli generate-s3-manifest --s3-base "${s3_base}"
    """
}
