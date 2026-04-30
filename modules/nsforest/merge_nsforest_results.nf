process merge_nsforest_results_process {
    tag "merge_${meta.organ}_${meta.first_author}_${meta.journal}_${meta.year}_${meta.embedding}_${meta.dataset_version_id}"
    label 'nsforest'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(meta),
          path(partial_csvs, stageAs: "partial_??.csv"),
          path(filtered_h5ad)

    output:
    tuple val(meta), path("*_results.csv"),                       emit: results_csv
    tuple val(meta), path("*_results_symbols.csv"),               emit: results_csv_symbols, optional: true
    tuple val(meta), path("*_results.pkl"),                       emit: results_pkl
    tuple val(meta), path("*_results_symbols.pkl"),               emit: results_pkl_symbols, optional: true
    tuple val(meta), path("*_markers.csv"),                       emit: markers
    tuple val(meta), path("*_markers_symbols.csv"),               emit: markers_symbols, optional: true
    tuple val(meta), path("*_markers_onTarget.csv"),              emit: markers_ontarget
    tuple val(meta), path("*_markers_onTarget_symbols.csv"),      emit: markers_ontarget_symbols, optional: true
    tuple val(meta), path("*_markers_onTarget_supp.csv"),         emit: markers_ontarget_supp
    tuple val(meta), path("*_markers_onTarget_supp_symbols.csv"), emit: markers_ontarget_supp_symbols, optional: true
    tuple val(meta), path("*_gene_selection.csv"),                emit: gene_selection
    tuple val(meta), path("*_gene_selection_symbols.csv"),        emit: gene_selection_symbols, optional: true

    script:
    """
    nsforest-cli merge-nsforest-results \
        --partial-files ${partial_csvs.join(',')} \
        --filtered-h5ad ${filtered_h5ad} \
        --cluster-header "${meta.author_cell_type}" \
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --journal "${meta.journal}" \
        --year "${meta.year}" \
        --embedding "${meta.embedding}" \
        --dataset-version-id "${meta.dataset_version_id}"
    """
}
