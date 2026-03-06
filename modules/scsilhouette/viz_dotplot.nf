process viz_dotplot_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'

    publishDir "${params.outdir}",
        mode: params.publish_mode,
        pattern: "outputs_*/**"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/*_dotplot_*.{html,svg}", optional: true),
          emit: plots

    script:
    """
    scsilhouette viz-dotplot \
        --h5ad-path ${h5ad} \
        --embedding-key "${meta.embedding}" \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}"
    """
}
