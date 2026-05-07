/**
 * Run NSForest Module (Parallelized by Cluster Batch)
 *
 * Runs NSForest algorithm to identify marker gene combinations.
 * Each cluster batch processed independently (one vs all).
 */
process run_nsforest_process {
    tag "run_nsforest_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}_${cluster}"
    label 'nsforest'

    input:
    tuple val(meta), path(h5ad), path(medians_csv), path(binary_scores_csv), val(cluster)

    output:
    tuple val(meta),
          path("*results*.csv"),
          emit: partial

    script:
    """
    # v2 - sanitized cluster filenames
    nsforest-cli run-nsforest \
        --h5ad-path ${h5ad} \
        --medians-csv ${medians_csv} \
        --binary-scores-csv ${binary_scores_csv} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
	--journal "${meta.journal}" \
        --year "${meta.year}" \
	--embedding "${meta.embedding}" \
	--dataset-version-id "${meta.dataset_version_id}" \
        --cluster-list '${cluster}' \
        --n-trees 1000 \
        --n-genes-eval 6
    """
}
