/**
 * Compare Published Module
 *
 * Compares our results with a published results CSV (for example the Jorstad et al. 2023 files on
 * Zenodo), cluster by cluster: markers, counts and metrics. Exact agreement is expected when the data
 * and the NS-Forest version are the same.
 *
 * Input:  tuple(meta, our results CSV, published CSV, label)
 * Output: comparison_<label>.csv (one row per cluster) and comparison_<label>.txt (the summary)
 */
process compare_published_process {
    tag "compare_${label}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(ours), path(published), val(label)

    output:
    tuple val(meta), path("comparison_*.{csv,txt}"), emit: comparison

    script:
    """
    nsforest-cli compare-published \
        --ours ${ours} \
        --published ${published} \
        --label "${label}"
    """
}
