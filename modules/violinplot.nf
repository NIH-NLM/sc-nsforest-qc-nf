process run_violinplot_process {

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
              path(base_sanitized_disease_tissue_medians_h5ad),
              path(base_sanitized_disease_tissue_binary_scores_h5ad),
	      path(base_sanitized_disease_tissue_nsforest_results_csv)
 
    output:
       tuple path(h5ad_file), val(label_key), val(embedding_key), val(organism), val(disease),
              val(filter), val(metric), val(save_scores), val(save_cluster_summary), val(save_annotation),
              val(tissue), val(author), val(publication_date), val(publication), val(cell_count),
              val(base),
              path(base_sanitized_h5ad),
              path(base_sanitized_disease_h5ad),
              path(base_sanitized_disease_tissue_h5ad),
              path(base_sanitized_disease_tissue_medians_h5ad),
              path(base_sanitized_disease_tissue_binary_scores_h5ad),
              path(base_sanitized_disease_tissue_nsforest_results_csv),
              path("${base}-sanitized-${disease}-${tissue}-violinplot.png"),
              path("${base}-sanitized-${disease}-${tissue}-violinplot.svg"),
              emit: violinplot_output_ch

    script:
    """
    nsforest-cli violinplot \
    --label-key $label_key \
    --png-out ${base}-sanitized-${disease}-${tissue}-violinplot.png \
    --svg-out ${base}-sanitized-${disease}-${tissue}-violinplot.svg \
    ${base_sanitized_disease_tissue_h5ad} \
    ${base_sanitized_disease_tissue_nsforest_results_csv} 
    """
}

