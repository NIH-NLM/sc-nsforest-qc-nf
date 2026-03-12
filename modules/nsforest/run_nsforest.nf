/**
 * Run NSForest Module (Parallelized by Cluster Batch)
 *
 * Runs NSForest algorithm to identify marker gene combinations.
 * Each cluster batch processed independently (one vs all).
 */
process run_nsforest_process {
    tag "run_nsforest_${meta.organ}_${meta.first_author}_${meta.year}_${cluster}"
    label 'nsforest'

    input:
    tuple val(meta), path(h5ad), path(medians_csv), path(binary_scores_csv), val(cluster)

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type.replace(' ','_')}_results_*.csv"),
          emit: partial

    script:
    """
    nsforest-cli run-nsforest \
        --h5ad-path ${h5ad} \
        --medians-csv ${medians_csv} \
        --binary-scores-csv ${binary_scores_csv} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}" \
        --cluster-list '${cluster}' \
        --n-trees 1000 \
        --n-genes-eval 6
    """
}
