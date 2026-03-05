process dendrogram_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode,
        pattern: "outputs_*/**"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path(h5ad),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/*.csv", optional: true),
          emit: stats
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/*.{csv,svg,html}", optional: true),
          emit: results

    script:
    """
    nsforest-cli dendrogram \
        --h5ad-path ${h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year}
    """
}
