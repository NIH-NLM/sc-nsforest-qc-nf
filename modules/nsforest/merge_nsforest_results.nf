process merge_nsforest_results_process {
    tag "merge_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta),
          path(partial_csvs, stageAs: "partial_??.csv"),
          path(filtered_h5ad)

    output:
    tuple val(meta),
          path("*results.csv"),
          path("*results.pkl"),
          path("*markers*.csv", optional: true),
          path("*gene_selection.csv", optional: true),
          emit: complete
    path "*results_symbols.csv", optional: true
    path "*results_symbols.pkl", optional: true
    path "*gene_selection_symbols.csv", optional: true

    script:
    """
    nsforest-cli merge-nsforest-results \
        --partial-files ${partial_csvs.join(',')} \
        --filtered-h5ad ${filtered_h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --journal "${meta.journal}" \
        --year "${meta.year}" \
        --embedding "${meta.embedding}" \
        --dataset-version-id "${meta.dataset_version_id}"
    """
}
