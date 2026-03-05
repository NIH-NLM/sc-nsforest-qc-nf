process viz_summary_process {
    tag "${meta.organ}_${meta.first_author}_${meta.year}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'

    publishDir "${params.outdir}",
        mode: params.publish_mode,
        pattern: "outputs_*/**"

    input:
    tuple val(meta),
          path(silhouette_scores),
          path(cluster_summary),
          path(annotation),
          path(nsforest_results)

    output:
    tuple val(meta),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_silhouette_fscore_summary.{html,svg}"),
          path("outputs_${meta.organ}_${meta.first_author}_${meta.year}/${meta.author_cell_type}_dataset_summary.csv"),
          emit: plots
      
    script:
    def fscore_flag       = nsforest_results.name != 'NO_FILE' ? "--fscore-path ${nsforest_results}" : ""
    def doi_flag          = meta.doi          ? "--doi \"${meta.doi}\""                   : ""
    def collection_flag   = meta.collection_name ? "--collection-name \"${meta.collection_name}\"" : ""
    def dataset_flag      = meta.dataset_title   ? "--dataset-title \"${meta.dataset_title}\""     : ""
    def journal_flag      = meta.journal         ? "--journal \"${meta.journal}\""                 : ""
    def coll_url_flag     = meta.collection_url  ? "--collection-url \"${meta.collection_url}\""   : ""
    def explorer_url_flag = meta.explorer_url    ? "--explorer-url \"${meta.explorer_url}\""       : ""
    def h5ad_url_flag     = meta.h5ad_url        ? "--h5ad-url \"${meta.h5ad_url}\""               : ""
    """
    scsilhouette viz-summary \
        --silhouette-score-path ${silhouette_scores} \
        --cluster-header "${meta.author_cell_type}" \
        --organ ${meta.organ} \
        --first-author ${meta.first_author} \
        --year ${meta.year} \
        --embedding-key ${meta.embedding} \
        ${doi_flag} \
        ${collection_flag} \
        ${dataset_flag} \
        ${journal_flag} \
        ${coll_url_flag} \
        ${explorer_url_flag} \
        ${h5ad_url_flag} \
        ${fscore_flag}
    """
}
