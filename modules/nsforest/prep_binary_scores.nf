process prep_binary_scores_process {
    tag "prep_binary_scores_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*binary_scores.csv"),         emit: csv
    tuple val(meta), path("*binary_scores_symbols.csv"), emit: csv_symbols, optional: true
    tuple val(meta), path("*binary_scores.pkl"),         emit: pkl
    tuple val(meta), path("*binary_scores_symbols.pkl"), emit: pkl_symbols, optional: true

    script:
    """
    nsforest-cli prep-binary-scores \
        --h5ad-path ${h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
	--journal "${meta.journal}" \
        --year "${meta.year}" \
	--embedding "${meta.embedding}" \
	--dataset-version-id "${meta.dataset_version_id}"
    """
}
