/**
 * Combinations Module
 *
 * For each cluster, evaluates the combinations of its own NS-Forest markers and its parent
 * subclass's markers with nsforest.evaluating.DecisionTree(combinations=True), and keeps the
 * combination with the best F-beta score (Jorstad et al. 2023, Zenodo 10.5281/zenodo.13172641).
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:            the cluster-level dataset (author_cell_type = the cluster header)
 *   - h5ad:            the h5ad both NS-Forest runs used
 *   - subclass_header: obs column of the parent subclasses
 *   - subclass_results: merged NSForest results CSV for the subclasses
 *   - cluster_results:  merged NSForest results CSV for the clusters
 *
 * Output:
 * -------
 * @emit results: tuple(meta, results_combinations_{prefix}.csv)
 */
process combinations_process {
    tag "combinations_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad), val(subclass_header), path(subclass_results), path(cluster_results)

    output:
    tuple val(meta), path("results_combinations_*.csv"), emit: results

    script:
    """
    nsforest-cli combinations \
        --h5ad-path ${h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --subclass-header "${subclass_header}" \
        --subclass-results ${subclass_results} \
        --cluster-results ${cluster_results} \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --journal "${meta.journal}" \
        --year "${meta.year}" \
        --embedding "${meta.embedding}" \
        --dataset-version-id "${meta.dataset_version_id}"
    """
}
