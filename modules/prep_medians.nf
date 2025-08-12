process prep_medians_process {

    tag "${h5ad_file}-${label_key}-${embedding_key}-${organism}-${disease}-${filter},${metric}-${save_scores}-${save_cluster_summary}-${save_annotation}-${tissue}-${author}-${publication_date}-${publication}-${cell_count}"

    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple path(h5ad_file), val(label_key), val(embedding_key), val(organism),val(disease),
	      val(filter), val(metric), val(save_scores), val(save_cluster_summary),val(save_annotation),
	      val(tissue), val(author), val(publication_date), val(publication),val(cell_count)

    output:
        tuple path(h5ad_file), val(label_key), val(embedding_key), val(organism),val(disease),
	      val(filter), val(metric), val(save_scores), val(save_cluster_summary),val(save_annotation),
	      val(tissue), val(author), val(publication_date), val(publication),val(cell_count),
	      path("medians*.h5ad"),
              emit: prep_medians_output_ch

    script:
    """
    nsforest-cli prep-medians --input-path $h5ad_file --cluster-header $label_key --output-path medians-$label_key-$embedding_key-$h5ad_file
    """
}

