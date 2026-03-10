process dendrogram_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path(h5ad),
          path("${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type.replace(' ','_')}_cluster_order.csv"),
          emit: stats
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type.replace(' ','_')}*.{csv,svg}"),
          emit: results

    script:
    """
    nsforest-cli dendrogram \
        --h5ad-path ${h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}"
    """
}
