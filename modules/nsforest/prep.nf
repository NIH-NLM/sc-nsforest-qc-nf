/* Prep: medians + binary scores in a single pass (one h5ad load / one densification) */
process prep_process {
    tag "prep_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}",
        mode: params.publish_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("medians_ensg*.csv"),          emit: medians_csv
    tuple val(meta), path("medians_symbols*.csv"),       emit: medians_csv_symbols, optional: true
    tuple val(meta), path("medians_ensg*.pkl"),          emit: medians_pkl
    tuple val(meta), path("medians_symbols*.pkl"),       emit: medians_pkl_symbols, optional: true
    tuple val(meta), path("binary_scores_ensg*.csv"),    emit: binary_csv
    tuple val(meta), path("binary_scores_symbols*.csv"), emit: binary_csv_symbols, optional: true
    tuple val(meta), path("binary_scores_ensg*.pkl"),    emit: binary_pkl
    tuple val(meta), path("binary_scores_symbols*.pkl"), emit: binary_pkl_symbols, optional: true

    script:
    """
    nsforest-cli prep \
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
