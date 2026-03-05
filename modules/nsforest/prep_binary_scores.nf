process prep_binary_scores_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'

    input:
    tuple val(meta), path(adata_prep)

    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/adata_prep.h5ad"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/*.{csv,svg,html,pkl}"),
          emit: complete

    script:
    """
    nsforest-cli prep-binary-scores \
        --h5ad-path ${adata_prep} \
        --cluster-header "${meta.author_cell_type}" \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year}
    """
}
