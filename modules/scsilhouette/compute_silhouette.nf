process compute_silhouette_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type.replace(' ','_')}*.{csv,json}"),
          emit: results

    script:
    """
    scsilhouette compute-silhouette \
        --h5ad-path ${h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --embedding-key "${meta.embedding}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}" \
        --disease "${meta.disease}" \
        --save-scores \
        --save-cluster-summary \
        --save-annotation
    """
}
