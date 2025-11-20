process sanitize_labels_process {

    tag "${h5ad_file.baseName}-${label_key}-${embedding_key}-${organism}-${disease}-${filter},${metric}-${save_scores}-${save_cluster_summary}-${save_annotation}-${tissue}-${author}-${publication_date}-${publication}-${cell_count}"

    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple path(h5ad_file), val(label_key), val(embedding_key), val(organism), val(disease),
              val(filter), val(metric), val(save_scores), val(save_cluster_summary), val(save_annotation),
              val(tissue), val(author), val(publication_date), val(publication), val(cell_count),
              val(base)

    output:
        tuple path(h5ad_file), val(label_key), val(embedding_key), val(organism), val(disease),
              val(filter), val(metric), val(save_scores), val(save_cluster_summary), val(save_annotation),
              val(tissue), val(author), val(publication_date), val(publication), val(cell_count),
              val(base),
              path("${base}-sanitized.h5ad"),
              emit: sanitize_output_ch

    script:
    """
    cp $h5ad_file ${base}-sanitized.h5ad
#    nsforest-cli sanitize-labels \
#    --h5ad-in=$h5ad_file \
#    --label-key=$label_key \
#    --h5ad-out=${base}-sanitized.h5ad
    """
}

