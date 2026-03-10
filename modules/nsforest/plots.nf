process plots_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad), path(results_csv)

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type.replace(' ','_')}*.{html,svg}"),
          emit: plots

    script:
    """
    nsforest-cli plots \
        --h5ad-path ${h5ad} \
        --results-csv ${results_csv} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}"
    """
}
