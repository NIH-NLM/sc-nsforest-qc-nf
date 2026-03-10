process prep_binary_scores_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'nsforest'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type.replace(' ','_')}_binary_scores.csv"),
          path("${meta.organ}_${meta.first_author}_${meta.year}_${meta.author_cell_type.replace(' ','_')}_binary_scores.pkl"),
          emit: complete

    script:
    """
    nsforest-cli prep-binary-scores \
        --h5ad-path ${h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}"
    """
}
