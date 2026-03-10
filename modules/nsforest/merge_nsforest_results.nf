process merge_nsforest_results_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(partial_csvs, stageAs: "partial_??.csv")

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type.replace(' ','_')}_results.csv"),
          path("${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type.replace(' ','_')}_results.pkl"),
          emit: complete

    script:
    """
    nsforest-cli merge-nsforest-results \
        --partial-files ${partial_csvs.join(',')} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}"
    """
}
