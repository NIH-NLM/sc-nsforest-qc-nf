process nsforest_process {

    tag "${h5ad_file}-${label_key}-${embedding_key}-${organism}-${disease}-${tissue}-${author}-${publication_date}-${publication}-${cell_count}"

    input:
       tuple path(h5ad_file), val(label_key), val(embedding_key), val(organism),
	      val(disease), val(tissue), val(author), val(publication_date), val(publication),
	      val(cell_count)

    output:
        tuple path(h5ad_file), val(label_key), val(embedding_key), val(organism), 
              val(disease), val(tissue), val(author), val(publication_date), val(publication),
	      val(cell_count), emit: nsforest_process_output_ch

    script:
    """
    nsforest.py --preprocess-adata-file -c $label_key $h5ad_file
    """
}

