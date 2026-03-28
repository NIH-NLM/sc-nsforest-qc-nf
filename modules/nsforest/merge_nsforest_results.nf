process merge_nsforest_results_process {
    tag "merge_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(partial_csvs, stageAs: "partial_??.csv")

    output:
    tuple val(meta),
          path("*results.csv"),
          path("*results.pkl"),
          path("*markers*.csv", optional: true),
          path("*gene_selection.csv", optional: true),
          emit: complete

    script:
    """
    nsforest-cli merge-nsforest-results \
        --partial-files ${partial_csvs.join(',')} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
	--journal "${meta.journal}" \
        --year "${meta.year}" \
        --embedding "${meta.embedding}" \
	--dataset-version-id "${meta.dataset_version_id}"
    """
}
