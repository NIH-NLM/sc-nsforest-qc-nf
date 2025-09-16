process prep_binary_scores_process {

    tag "${h5ad_file}-${label_key}-${embedding_key}-${organism}-${disease}-${filter},${metric}-${save_scores}-${save_cluster_summary}-${save_annotation}-${tissue}-${author}-${publication_date}-${publication}-${cell_count}"

    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple path(h5ad_file), val(label_key), val(embedding_key), val(organism), val(disease),
              val(filter), val(metric), val(save_scores), val(save_cluster_summary), val(save_annotation),
              val(tissue), val(author), val(publication_date), val(publication), val(cell_count),
              val(base),
              path(base_sanitized_h5ad),
              path(base_sanitized_disease_h5ad),
              path(base_sanitized_disease_tissue_h5ad),
              path(base_sanitized_disease_tissue_medians_h5ad)

     output:
        tuple path(h5ad_file), val(label_key), val(embedding_key), val(organism), val(disease),
              val(filter), val(metric), val(save_scores), val(save_cluster_summary), val(save_annotation),
              val(tissue), val(author), val(publication_date), val(publication), val(cell_count),
              val(base),
              path(base_sanitized_h5ad),
              path(base_sanitized_disease_h5ad),
              path(base_sanitized_disease_tissue_h5ad),
              path(base_sanitized_disease_tissue_medians_h5ad),
              path("${base}-sanitized-${disease}-${tissue}-binary-scores.h5ad"),
              emit: prep_binary_scores_output_ch

    script:
    """
    nsforest-cli prep-binary-scores \
    --h5ad-in=$base_sanitized_disease_tissue_medians_h5ad \
    --label-key=$label_key \
    --h5ad-out=${base}-sanitized-${disease}-${tissue}-binary-scores.h5ad
    """
}

