/**
 * Viz Summary Module
 *
 * Generates an interactive silhouette F-score summary plot combining
 * silhouette scores with NSForest F-scores per cluster.
 * Also writes a dataset-level summary CSV.
 *
 * Input:
 * ------
 * @param tuple:
 *   - meta:             Map with organ, first_author, year, author_cell_type, embedding, doi, etc.
 *   - silhouette_scores:{prefix}_silhouette_scores.csv
 *   - cluster_summary:  {prefix}_cluster_summary.csv
 *   - annotation:       {prefix}_annotation.json
 *   - nsforest_results: {prefix}_results.csv (or NO_FILE sentinel)
 *
 * Output:
 * -------
 * @emit plots: tuple(meta, [summary SVG/HTML/CSV, dataset_summary CSV])
 *   Flat filenames: {organ}_{first_author}_{year}_{cluster_header_safe}_*.{csv,svg,html,json}
 */
process viz_summary_process {
    tag "viz_summary_${meta.organ}_${meta.first_author}_${meta.year}"
    label 'scsilhouette'
    containerOptions '--entrypoint ""'
    publishDir "${params.outdir}",
        mode: params.publish_mode,
        pattern: "*.{csv,svg,html,json}"

    input:
    tuple val(meta),
          path(silhouette_scores),
          path(cluster_summary),
          path(annotation),
          path(nsforest_results)

    output:
    tuple val(meta),
          path("${meta.organ}_${meta.first_author}_${meta.year}_*.{csv,svg,html,json,log}", optional: true),
          emit: plots

    script:
    def fscore_flag       = nsforest_results.name != 'NO_FILE' ? "--fscore-path ${nsforest_results}" : ""
    def doi_flag          = meta.doi             ? "--doi \"${meta.doi}\""                         : ""
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
        --organ "${meta.organ}" \
        --first-author "${meta.first_author}" \
        --year "${meta.year}" \
        --embedding-key "${meta.embedding}" \
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
