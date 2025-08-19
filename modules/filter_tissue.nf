process filter_tissue_process {

    tag "${h5ad_file.baseName}-${label_key}-${embedding_key}-${organism}-${disease}-${filter},${metric}-${save_scores}-${save_cluster_summary}-${save_annotation}-${tissue}-${author}-${publication_date}-${publication}-${cell_count}"

    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple path(h5ad_file), val(label_key), val(embedding_key), val(organism), val(disease),
              val(filter), val(metric), val(save_scores), val(save_cluster_summary), val(save_annotation),
              val(tissue), val(author), val(publication_date), val(publication), val(cell_count),
              val(base),
              path("${base}-sanitized.h5ad,
              path("${base}-sanitized-${disease}.h5ad")
    output:
        tuple path(h5ad_file), val(label_key), val(embedding_key), val(organism), val(disease),
              val(filter), val(metric), val(save_scores), val(save_cluster_summary), val(save_annotation),
              val(tissue), val(author), val(publication_date), val(publication), val(cell_count),
              val(base),
              path("${base}-sanitized.h5ad,
              path("${base}-sanitized-${disease}.h5ad"),
              path("${base}-sanitized-${disease}-${tissue}.h5ad"),
              emit: filter_kidney_output_ch

    script:
    """
    nsforest-cli filter-by-obs \
    --obs-key tissue \
    --value ${tissue} \
    ${base}-sanitized-${disease}.h5ad \
    ${base}-sanitized-${disease}-${tissue}.h5ad
    """
}

